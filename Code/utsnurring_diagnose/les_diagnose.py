"""Leser utsnurring_diagnose over USB: begge flankene fra Hall-latchen, merket S og N.

    python3 les_diagnose.py                  # finner porten selv
    python3 les_diagnose.py --stopp 5

Skriver hver flanke med turtallet regnet fra forrige flanke med samme pol, og markerer
med <<< når to like poler kommer etter hverandre (sensoren gikk glipp av den andre magneten).
Statuslinjene fra ESP-en vises bare når det ikke kommer flanker, så du ser hvilket nivå
utgangen står fast på. Lagres i data/diagnose_<tid>.csv (t_s,pol). Krever pyserial.
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


def oppsummer(rader):
    s = [t for t, p in rader if p == "S"]
    n = [t for t, p in rader if p == "N"]
    dobbel = [(i, p) for i, ((_, p), (_, q)) in enumerate(zip(rader, rader[1:]), 1) if p == q]
    print(f"\n{len(rader)} flanker: S = {len(s)}, N = {len(n)}, like etter hverandre = {len(dobbel)}")
    for navn, tt in (("S", s), ("N", n)):
        if len(tt) >= 2:
            print(f"  {navn}: siste flanke t = {tt[-1]:.4f} s, turtall der {60 / (tt[-1] - tt[-2]):.0f} o/min")
    if rader:
        siste = rader[-1][1]
        print(f"  Utgangen sto til slutt {'lav (sydpol sist sett)' if siste == 'S' else 'høy (nordpol sist sett)'}.")
    if dobbel:
        print("  Første gang to like kom etter hverandre: flanke", dobbel[0][0] + 1, dobbel[0][1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port")
    ap.add_argument("--stopp", type=float, default=5.0, help="sekunder uten flanker før målingen avsluttes")
    a = ap.parse_args()
    port = a.port or finn_port()

    rader, sist_t = [], {}
    sist = None
    with serial.Serial(port, 115200, timeout=0.2) as s:
        s.reset_input_buffer()
        print(f"Lytter på {port}. Snurr skiva (Ctrl+C for å avslutte).")
        try:
            while True:
                linje = s.readline().decode(errors="replace").strip()
                if linje.startswith("#"):
                    if sist is None or time.monotonic() - sist > 1.0:
                        print("  ", linje)
                elif "," in linje:
                    try:
                        ts, pol = linje.split(",")
                        t = float(ts)
                    except ValueError:
                        continue  # overskriften
                    if rader and t - rader[-1][0] > a.stopp:
                        rader, sist_t = [], {}  # gammel flanke fra før porten ble åpnet
                    merke = "  <<< samme pol som forrige" if rader and rader[-1][1] == pol else ""
                    rpm = f"{60 / (t - sist_t[pol]):7.0f} o/min" if pol in sist_t else "      -     "
                    sist_t[pol] = t
                    rader.append((t, pol))
                    sist = time.monotonic()
                    print(f"  {len(rader):5d}  {t:11.4f}  {pol}  {rpm}{merke}")
                if rader and sist is not None and time.monotonic() - sist > a.stopp:
                    break
        except KeyboardInterrupt:
            pass

    oppsummer(rader)
    if rader:
        DATA.mkdir(exist_ok=True)
        fil = DATA / f"diagnose_{datetime.now():%Y%m%d_%H%M%S}.csv"
        fil.write_text("t_s,pol\n" + "".join(f"{t:.6f},{p}\n" for t, p in rader))
        print("Lagret", fil)


if __name__ == "__main__":
    main()
