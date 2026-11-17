"""Leser utsnurring_diagnose over USB eller BLE: begge flankene fra Hall-latchen, merket S og N.

    python3 les_diagnose.py                  # USB, finner porten selv
    python3 les_diagnose.py --ble            # trådløst, kobler til "KERS-utsnurring"
    python3 les_diagnose.py --ble --stopp 600
    python3 les_diagnose.py --ble --plott    # med graf av turtallet som tegnes live

Skriver hver flanke med turtallet regnet fra forrige flanke med samme pol, og markerer
med <<< når to like poler kommer etter hverandre (sensoren gikk glipp av den andre magneten).
Statuslinjene fra ESP-en vises bare når det ikke kommer flanker, så du ser hvilket nivå
utgangen står fast på. Lagres i data/diagnose_<tid>.csv (t_s,pol).
USB krever pyserial, BLE krever bleak.
"""

import argparse
import asyncio
import glob
import threading
import time
from datetime import datetime
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"

BLE_NAVN = "KERS-utsnurring"
FLANKER_UUID = "4a4bf346-512e-49b7-97a3-f5c3cd2b8a71"
STATUS_UUID = "c0ebe3b8-889d-414e-87a1-23cde9e22c61"


class Logg:
    """Samler flankene og skriver dem fortløpende, uansett kilde."""

    def __init__(self, stopp):
        self.stopp = stopp
        self.rader, self.sist_t = [], {}
        self.punkter = {"S": [], "N": []}  # (t, o/min) per pol, til grafen
        self.sist = None  # time.monotonic() ved siste flanke
        self.siste_status = ""
        self.avbrutt = threading.Event()  # settes når grafvinduet lukkes

    def flanke(self, t, pol):
        if self.rader and t - self.rader[-1][0] > self.stopp:
            self.rader, self.sist_t = [], {}  # gammel flanke fra før målingen startet
            self.punkter = {"S": [], "N": []}
        merke = "  <<< samme pol som forrige" if self.rader and self.rader[-1][1] == pol else ""
        rpm = f"{60 / (t - self.sist_t[pol]):7.0f} o/min" if pol in self.sist_t else "      -     "
        if pol in self.sist_t and t > self.sist_t[pol]:
            self.punkter[pol].append((t, 60 / (t - self.sist_t[pol])))
        self.sist_t[pol] = t
        self.rader.append((t, pol))
        self.sist = time.monotonic()
        print(f"  {len(self.rader):5d}  {t:11.4f}  {pol}  {rpm}{merke}")

    def status(self, linje):
        self.siste_status = linje
        if self.sist is None or time.monotonic() - self.sist > 1.0:
            print("  ", linje)

    def ferdig(self):
        return self.avbrutt.is_set() or bool(self.rader) and self.sist is not None and time.monotonic() - self.sist > self.stopp


def finn_port():
    porter = sorted(glob.glob("/dev/cu.usbmodem*") + glob.glob("/dev/cu.usbserial*") + glob.glob("/dev/ttyACM*"))
    if not porter:
        raise SystemExit("Fant ingen USB-port. Er ESP-en koblet til?")
    return porter[0]


def les_usb(logg, port):
    import serial

    with serial.Serial(port, 115200, timeout=0.2) as s:
        s.reset_input_buffer()
        print(f"Lytter på {port}. Snurr skiva (Ctrl+C for å avslutte).")
        while not logg.ferdig():
            linje = s.readline().decode(errors="replace").strip()
            if linje.startswith("#"):
                logg.status(linje)
            elif "," in linje:
                try:
                    ts, pol = linje.split(",")
                    logg.flanke(float(ts), pol)
                except ValueError:
                    continue  # overskriften


async def les_ble(logg, navn):
    from bleak import BleakClient, BleakScanner

    print(f"Leter etter {navn} ...")
    enhet = await BleakScanner.find_device_by_name(navn, timeout=15)
    if enhet is None:
        raise SystemExit(f"Fant ikke {navn}. Står ESP-en på, og er Bluetooth slått på?")

    forrige, omslag = None, 0  # micros() går rundt etter 2^32 µs, omtrent 71 min

    def flanker(_, data):
        nonlocal forrige, omslag
        for k in range(0, len(data) - 4, 5):
            us = int.from_bytes(data[k:k + 4], "little")
            if forrige is not None and us < forrige - 2**31:
                omslag += 1
            forrige = us
            logg.flanke((us + omslag * 2**32) / 1e6, chr(data[k + 4]))

    def status(_, data):
        logg.status(data.decode(errors="replace"))

    while not logg.ferdig():
        frakoblet = asyncio.Event()
        async with BleakClient(enhet, disconnected_callback=lambda _: frakoblet.set()) as klient:
            await klient.start_notify(FLANKER_UUID, flanker)
            await klient.start_notify(STATUS_UUID, status)
            print(f"Koblet til {navn} ({enhet.address}). Snurr skiva (Ctrl+C for å avslutte).")
            while not logg.ferdig() and not frakoblet.is_set():
                await asyncio.sleep(0.1)
        if not logg.ferdig():
            print("   Mistet forbindelsen, kobler til igjen ...")
            await asyncio.sleep(0.5)


def vis_graf(logg, kilde):
    """Kjører kilde() i en egen tråd og tegner turtallet mot tid mens dataene kommer inn.

    Når målingen er ferdig, blir vinduet stående. Lukk det for å lagre.
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    feil = []

    def kjor():
        try:
            kilde()
        except BaseException as e:  # SystemExit fra «fant ikke enheten» o.l.
            feil.append(e)

    traad = threading.Thread(target=kjor, daemon=True)
    traad.start()

    fig, ax = plt.subplots(figsize=(9, 5))
    fig.canvas.manager.set_window_title("Utsnurring live")
    linjer = {p: ax.plot([], [], ".", ms=3, label=f"{p}-pol")[0] for p in "SN"}
    ax.set_xlabel("Tid fra første flanke [s]")
    ax.set_ylabel("Turtall [o/min]")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")

    def oppdater(_):
        if logg.rader:
            t0 = logg.rader[0][0]
            t_maks, w_maks = 1.0, 100.0
            for p, linje in linjer.items():
                pk = list(logg.punkter[p])
                linje.set_data([t - t0 for t, _ in pk], [w for _, w in pk])
                if pk:
                    t_maks = max(t_maks, pk[-1][0] - t0)
                    w_maks = max(w_maks, max(w for _, w in pk))
            # Settes for hånd: set_ylim(bottom=0) skrur av autoskaleringen, og da ble grafen stående tom.
            ax.set_xlim(0, t_maks * 1.05)
            ax.set_ylim(0, w_maks * 1.1)
        alle = logg.punkter["S"][-1:] + logg.punkter["N"][-1:]
        na = max(alle)[1] if alle else 0
        tekst = f"{na:.0f} o/min   {len(logg.rader)} flanker"
        tapt = logg.siste_status.partition("tapt=")[2]
        if tapt and tapt != "0":
            tekst += f"   ({tapt} flanker tapt)"
        if not traad.is_alive():
            tekst += "   FERDIG, lukk vinduet for å lagre" if not feil else f"   STOPPET: {feil[0]}"
        ax.set_title(tekst)

    animasjon = FuncAnimation(fig, oppdater, interval=200, cache_frame_data=False)  # noqa: F841
    try:
        plt.show()
    except KeyboardInterrupt:
        pass
    logg.avbrutt.set()
    traad.join(timeout=5)
    if feil and not isinstance(feil[0], KeyboardInterrupt):
        print(feil[0])


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
    ap.add_argument("--ble", action="store_true", help="les trådløst over BLE i stedet for USB")
    ap.add_argument("--navn", default=BLE_NAVN, help="BLE-navnet til ESP-en")
    ap.add_argument("--plott", action="store_true", help="vis turtallet i en graf som oppdateres live")
    ap.add_argument("--stopp", type=float, default=5.0, help="sekunder uten flanker før målingen avsluttes")
    a = ap.parse_args()

    logg = Logg(a.stopp)
    if a.ble:
        kilde = lambda: asyncio.run(les_ble(logg, a.navn))
    else:
        port = a.port or finn_port()
        kilde = lambda: les_usb(logg, port)
    try:
        if a.plott:
            vis_graf(logg, kilde)
        else:
            kilde()
    except KeyboardInterrupt:
        pass

    oppsummer(logg.rader)
    if logg.rader:
        DATA.mkdir(exist_ok=True)
        fil = DATA / f"diagnose_{datetime.now():%Y%m%d_%H%M%S}.csv"
        fil.write_text("t_s,pol\n" + "".join(f"{t:.6f},{p}\n" for t, p in logg.rader))
        print("Lagret", fil)


if __name__ == "__main__":
    main()
