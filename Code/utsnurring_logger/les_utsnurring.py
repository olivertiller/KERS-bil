"""Leser pulstider fra utsnurring_logger (ESP32/Arduino) over USB og lagrer dem som CSV.

    python3 les_utsnurring.py                 # finner porten selv
    python3 les_utsnurring.py --port /dev/cu.usbmodem1101 --stopp 5

Viser turtallet for hver omdreining mens svinghjulet spinner. Når det har kommet pulser og
det så har gått --stopp sekunder uten nye, lagres målingen i data/utsnurring_<tid>.csv
(én kolonne t_s, én rad per omdreining), klar for  python3 ../kers_sim.py --utsnurring fil.csv
Ctrl+C lagrer også det som er målt så langt.

Krever pyserial:  python3 -m pip install pyserial
"""

import argparse
import glob
import time
from datetime import datetime
from pathlib import Path

import serial

DATA = Path(__file__).resolve().parent / "data"


def finn_port():
    porter = sorted(glob.glob("/dev/cu.usbmodem*") + glob.glob("/dev/cu.usbserial*") + glob.glob("/dev/ttyACM*"))
    if not porter:
        raise SystemExit("Fant ingen USB-port. Er ESP-en koblet til?")
    return porter[0]


def lagre(tider):
    DATA.mkdir(exist_ok=True)
    fil = DATA / f"utsnurring_{datetime.now():%Y%m%d_%H%M%S}.csv"
    fil.write_text("t_s\n" + "".join(f"{t:.6f}\n" for t in tider))
    return fil


def oppsummer(tider):
    if len(tider) < 3:
        print(f"\nBare {len(tider)} pulser. Sjekk at én magnet med sydpol og én med nordpol "
              "passerer skriftsiden på sensoren (en latch trenger begge for å gi pulser).")
        return
    rpm = [60.0 / (b - a) for a, b in zip(tider, tider[1:])]
    print(f"\n{len(tider)} pulser over {tider[-1] - tider[0]:.1f} s, maks {max(rpm):.0f} o/min, "
          f"siste {rpm[-1]:.0f} o/min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", help="seriellport (standard: første /dev/cu.usbmodem*)")
    ap.add_argument("--stopp", type=float, default=5.0, help="sekunder uten puls før målingen avsluttes")
    a = ap.parse_args()
    port = a.port or finn_port()

    tider = []
    sist = None
    with serial.Serial(port, 115200, timeout=0.2) as s:
        s.reset_input_buffer()
        print(f"Lytter på {port}. Spinn opp svinghjulet (Ctrl+C for å avslutte).")
        try:
            while True:
                linje = s.readline().decode(errors="replace").strip()
                if linje:
                    try:
                        t = float(linje)
                    except ValueError:
                        continue  # overskriften "t_s" og annen tekst
                    if tider and t <= tider[-1]:
                        print("Tiden gikk bakover (ESP-en startet på nytt?), begynner på nytt.")
                        tider = []
                    elif tider and t - tider[-1] > a.stopp:
                        # Gammel puls som lå i ESP-ens buffer da porten ble åpnet, eller en
                        # enkeltpuls lenge før spinnet: start serien på nytt.
                        tider = []
                    tider.append(t)
                    sist = time.monotonic()
                    if len(tider) >= 2:
                        print(f"  {len(tider):5d}  {60.0 / (tider[-1] - tider[-2]):7.0f} o/min")
                    else:
                        print("  første puls")
                elif tider and time.monotonic() - sist > a.stopp:
                    break
        except KeyboardInterrupt:
            pass

    oppsummer(tider)
    if tider:
        print("Lagret", lagre(tider))


if __name__ == "__main__":
    main()
