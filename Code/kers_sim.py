"""
Simulering av KERS-bilen: svinghjul -> girutveksling -> drivhjul.

Finner utvekslingen i = omega_svinghjul / omega_hjul som gir lengst kjøredistanse,
og kan tilpasse tapsmodellen for svinghjulet til en målt utsnurringskurve.

Bruk:
    python3 kers_sim.py                          # sveip over i med antatte tapstall
    python3 kers_sim.py --utsnurring maling.csv  # tilpass tapene til målt utsnurring først
    python3 kers_sim.py --rpm 3500 --frihjul     # annet startturtall / med frihjul

Utsnurringsfila er en CSV med overskriftsrad og enten
    t_s,pol:      diagnosefil fra Code/utsnurring_diagnose (S- og N-flanker, to per omdreining)
    én kolonne:   tidspunkt [s] for hver puls fra Code/utsnurring_logger (én per omdreining)
    to kolonner:  t [s], turtall [o/min]
Flankefilene kan ha mange trekk etter hverandre: hvert trekk kuttes ved toppen, og alle tilpasses\nsamtidig. Flere filer kan gis etter --utsnurring.

Modell:
    Svinghjulets tapsmoment (lager + luft):  T_tap = c0 + c1*w + c2*w^2
    Rullemotstand:                           F_rr  = C_rr*m*g
    Utvekslingen har virkningsgrad eta per trinn.
    Bilen holdes løftet under opptrekk og settes ned med svinghjulet i gang. Drivhjulene
    slurer da til farten til bil og svinghjul er lik, deretter ruller de uten slipp.
"""

import argparse
import math
from dataclasses import dataclass, replace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

g = 9.81


@dataclass
class Bil:
    I: float = 6.0e-5 + 6.5e-8  # svinghjul + aksel [kg m^2]
    m: float = 0.25  # bilens totale masse [kg] (antatt)
    r_h: float = 0.030  # drivhjulradius [m]
    J_hjul: float = 0.0  # samlet treghetsmoment for alle hjul [kg m^2]
    C_rr: float = 0.03  # rullemotstandskoeffisient (antatt)
    eta_trinn: float = 0.95  # virkningsgrad per tannhjulstrinn (antatt)
    n_trinn: int = 2  # antall tannhjulstrinn
    mu: float = 0.6  # friksjon mellom dekk og underlag (antatt)
    andel_drivhjul: float = 0.5  # andel av vekten som hviler på drivhjulene
    c0: float = 2e-4  # konstant friksjonsmoment i lagrene [N m] (antatt, mål dette!)
    c1: float = 0.0  # viskøst tap [N m s/rad]
    c2: float = 3.5e-10  # luftmotstand på skiva [N m s^2/rad^2] (overslag)
    frihjul: bool = False  # frihjul mellom svinghjul og drivhjul


def tapsmoment(bil, w):
    return bil.c0 + bil.c1 * w + bil.c2 * w * w if w > 0 else 0.0


def simuler(bil, i, rpm0, dt=1e-3, t_maks=900.0):
    """Simulerer én kjøring. Returnerer tidsserier og nøkkeltall."""
    eta = bil.eta_trinn**bil.n_trinn
    k = bil.r_h / i  # v = k*w når drivhjulene ruller uten slipp
    m_eff = bil.m + bil.J_hjul / bil.r_h**2
    F_rr = bil.C_rr * bil.m * g
    F_maks = bil.mu * bil.andel_drivhjul * bil.m * g
    tol = 1e-5

    def koblet(T_tap, F_rr_t, kk):
        # Ruller uten slipp: felles akselerasjon. kk = 1/eta når effekten går fra
        # svinghjul til hjul, kk = eta når bilen driver svinghjulet.
        w_dot = -(T_tap + F_rr_t * k * kk) / (bil.I + m_eff * k * k * kk)
        return w_dot, m_eff * k * w_dot + F_rr_t

    w, v, s, t = rpm0 * 2 * math.pi / 60, 0.0, 0.0, 0.0
    E0 = 0.5 * bil.I * w**2
    logg = [(t, w, v, s)]
    n = 0
    while t < t_maks:
        T_tap = tapsmoment(bil, w)
        slipp = k * w - v  # drivhjulenes periferihastighet minus bilens fart
        F_rr_t = F_rr if v > 0 else 0.0

        # F_d er drivkraften fra drivhjulene på bilen
        if abs(slipp) <= tol:
            w_dot, F_d = koblet(T_tap, F_rr_t, 1 / eta)
            if F_d < 0:
                if bil.frihjul:  # frihjulet løser ut, svinghjul og bil går hver for seg
                    w_dot, F_d = -T_tap / bil.I, 0.0
                else:
                    w_dot, F_d = koblet(T_tap, F_rr_t, eta)
            ruller = abs(F_d) <= F_maks
            if not ruller:  # for stor kraft, hjulene begynner å spinne
                F_d = math.copysign(F_maks, F_d)
        else:
            ruller = False
            F_d = 0.0 if (slipp < 0 and bil.frihjul) else math.copysign(F_maks, slipp)

        if not ruller:
            T_g = F_d * k / eta if F_d > 0 else F_d * k * eta
            w_dot = -(T_tap + T_g) / bil.I
        v_dot = (F_d - F_rr_t) / m_eff

        w_ny = max(w + w_dot * dt, 0.0)
        v_ny = max(v + v_dot * dt, 0.0)
        # Slipp som skifter fortegn betyr at hjul og bil har nådd samme fart:
        # lås dem sammen med bevart (redusert) bevegelsesmengde, som et uelastisk støt.
        if (k * w - v) * (k * w_ny - v_ny) < 0 and not (bil.frihjul and k * w_ny < v_ny):
            v_ny = (m_eff * v_ny + bil.I * w_ny / k) / (m_eff + bil.I / k**2)
            w_ny = v_ny / k
        w, v = w_ny, v_ny
        s += v * dt
        t += dt
        n += 1
        if n % 20 == 0:
            logg.append((t, w, v, s))
        if v <= 0 and k * w <= tol:
            break

    logg.append((t, w, v, s))
    t_a, w_a, v_a, s_a = map(np.array, zip(*logg))
    return {"t": t_a, "w": w_a, "v": v_a, "s": s_a, "distanse": s, "tid": t,
            "v_maks": v_a.max(), "E0": E0}


def les_utsnurring(fil):
    """Leser en utsnurringsfil og returnerer en liste med snurr [(t, turtall), ...], hvert fra toppen."""
    with open(fil) as f:
        hode = f.readline().strip().lower().split(",")
    if "pol" in hode:  # diagnosefil fra utsnurring_diagnose: t_s,pol (S og N, to flanker per omdreining)
        t = np.loadtxt(fil, delimiter=",", skiprows=1, usecols=0, ndmin=1)
        return utsnurring_fra_flanker(t, 2)
    data = np.loadtxt(fil, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[1] == 1:  # pulstider fra utsnurring_logger (én puls per omdreining)
        return utsnurring_fra_flanker(data[:, 0], 1)
    return [(data[:, 0], data[:, 1])]


KUTT_ETTER_TOPP = 0.3  # s; de første tiendelene etter toppen bremser snora/løkka ekstra (målt 2026-11-06)


def utsnurring_fra_flanker(t, per_omdr, min_rpm=300.0, min_punkter=10, kutt=KUTT_ETTER_TOPP):
    """Snurr fra flanketider. Hvert punkt bruker én hel omdreining (flanke k til k + per_omdr), så små
    feil i magnetplasseringen ikke gir rippel. Fila deles der det er pause mellom flankene (flere trekk
    i samme fil). Hvert snurr starter `kutt` s etter toppen, der drillen slapp eller snora gikk av
    (overgangen rett etter at løkka slipper bremser ekstra og er ikke en del av tapsmodellen). Snurr som
    aldri når min_rpm, eller der turtallet stiger igjen underveis (dreiing for hånd), hoppes over."""
    pause = 2.0 / per_omdr  # s; tilsvarer under 30 o/min
    serier = []
    for b in np.split(t, np.where(np.diff(t) > pause)[0] + 1):
        if len(b) <= per_omdr + min_punkter:
            continue
        dt = b[per_omdr:] - b[:-per_omdr]
        tm, rpm = 0.5 * (b[per_omdr:] + b[:-per_omdr]), 60.0 / dt
        k = int(np.argmax(rpm))
        k = int(np.searchsorted(tm, tm[k] + kutt))
        r = rpm[k:]
        frisnurr = np.mean(np.diff(r) > 0.005 * r[:-1]) < 0.15  # synker jevnt, ingen hånd på hjulet
        if rpm[k] >= min_rpm and len(r) >= min_punkter and frisnurr:
            serier.append((tm[k:], r))
    if not serier:
        raise SystemExit(f"Fant ingen utsnurring over {min_rpm:.0f} o/min i fila.")
    return serier


def _modell(p, t, w0, I):
    f = lambda _, y: [-(p[0] + p[1] * y[0] + p[2] * y[0] ** 2) / I if y[0] > 0 else 0.0]
    return solve_ivp(f, (0, t[-1]), [w0], t_eval=t, rtol=1e-8, atol=1e-6).y[0]


def tilpass_utsnurring(serier, I):
    """Tilpasser c0, c1, c2 slik at I*dw/dt = -(c0 + c1*w + c2*w^2) treffer alle snurrene samtidig.
    Hvert snurr starter fra sitt eget toppturtall. Returnerer (c0, c1, c2) og modellkurvene."""
    data = [(t - t[0], rpm * 2 * math.pi / 60) for t, rpm in serier]
    skala = np.array([1e-3, 1e-6, 1e-8])  # størrelsesorden for målt tap; for små skalaer satt fast
    rest = lambda q: np.concatenate([_modell(q * skala, t, w[0], I) - w for t, w in data])
    # Startgjetning fra snittmomentet i målingene. En fast gjetning på ~1 mNm stopper modellskiva
    # altfor tidlig når tapet er lite (SKF-lagrene, 0,2 mNm), og da blir løseren sittende fast.
    t_snitt = np.mean([I * (w[0] - w[-1]) / t[-1] for t, w in data])
    x0 = [0.7 * t_snitt / skala[0], 0.15 * t_snitt / 100 / skala[1], 0.15 * t_snitt / 100**2 / skala[2]]
    res = least_squares(rest, x0=x0, bounds=(0, np.inf))
    p = res.x * skala
    return p, [_modell(p, t, w[0], I) for t, w in data]


def tid_mellom(t, rpm, hoy, lav):
    """Tiden snurret bruker fra turtallet hoy til lav (lineær interpolasjon), eller nan."""
    if rpm[0] < hoy or rpm[-1] > lav:
        return float("nan")
    t_h = np.interp(-hoy, -rpm, t)  # rpm synker, så snu fortegnet for interp
    t_l = np.interp(-lav, -rpm, t)
    return t_l - t_h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rpm", type=float, default=4000, help="svinghjulets turtall etter opptrekk")
    ap.add_argument("--utsnurring", nargs="+", help="CSV-fil(er) med målt utsnurring; alle trekk tilpasses samtidig")
    ap.add_argument("--utelat", type=int, nargs="+", default=[], help="trekknumre (fra tabellen) som ikke skal brukes")
    ap.add_argument("--frihjul", action="store_true", help="simuler med frihjul")
    ap.add_argument("--i-min", type=float, default=2)
    ap.add_argument("--i-maks", type=float, default=40)
    ap.add_argument("--figur", default="kers_sim.png")
    a = ap.parse_args()

    bil = Bil(frihjul=a.frihjul)

    if a.utsnurring:
        serier = [sn for fil in a.utsnurring for sn in les_utsnurring(fil)]
        nummer = [n for n in range(1, len(serier) + 1) if n not in a.utelat]
        if a.utelat:
            print(f"Utelater trekk {', '.join(map(str, a.utelat))}")
            serier = [serier[n - 1] for n in nummer]
        (c0, c1, c2), w_fit = tilpass_utsnurring(serier, bil.I)
        bil = replace(bil, c0=c0, c1=c1, c2=c2)
        hoy, lav = 1000, 300  # felles turtallsvindu for å sammenligne trekkene
        print(f"{len(serier)} utsnurringer fra {len(a.utsnurring)} fil(er):")
        print(f"  {'nr':>3} {'topp [o/min]':>13} {'til stopp [s]':>14} {f'{hoy}->{lav} [s]':>13} {'snittmoment [mNm]':>18}")
        tider = []
        for n, (t, rpm) in zip(nummer, serier):
            dt = tid_mellom(t, rpm, hoy, lav)
            tider.append(dt)
            moment = bil.I * (hoy - lav) * 2 * math.pi / 60 / dt * 1e3
            print(f"  {n:3d} {rpm[0]:13.0f} {t[-1] - t[0]:14.2f} {dt:13.3f} {moment:18.2f}")
        tider = np.array(tider)
        ok = tider[np.isfinite(tider)]
        if len(ok) >= 2:
            snitt, sd = ok.mean(), ok.std(ddof=1)
            m = bil.I * (hoy - lav) * 2 * math.pi / 60 / snitt * 1e3
            print(f"  Snitt {hoy}->{lav} o/min: {snitt:.3f} s ± {sd:.3f} s (standardavvik, {sd / snitt * 100:.1f} %), "
                  f"snittmoment {m:.2f} mNm")
        print(f"Tilpasset tapsmodell (alle samtidig): c0 = {c0:.3e} N m, c1 = {c1:.3e} N m s, c2 = {c2:.3e} N m s^2")
        topp = max(rpm[0] for _, rpm in serier)
        if a.rpm > 1.1 * topp:
            print(f"OBS: --rpm {a.rpm:.0f} er over høyeste målte turtall ({topp:.0f}). Tapsmodellen ekstrapoleres.")
        fig, ax = plt.subplots()
        for n, ((t, rpm), w) in enumerate(zip(serier, w_fit)):
            l, = ax.plot(t - t[0], rpm, ".", ms=2, alpha=0.6)
            ax.plot(t - t[0], w * 60 / (2 * math.pi), color=l.get_color(), lw=0.8, label="modell" if n == 0 else None)
        ax.set(xlabel="Tid fra toppen [s]", ylabel="Turtall [o/min]",
               title=f"Utsnurring av svinghjul, {len(serier)} trekk (punkter = målt)")
        ax.legend()
        fig.savefig("utsnurring_tilpasning.png", dpi=150, bbox_inches="tight")

    utvekslinger = np.arange(a.i_min, a.i_maks + 1e-9, 1.0)
    resultater = [simuler(bil, i, a.rpm) for i in utvekslinger]
    distanser = np.array([r["distanse"] for r in resultater])
    beste = int(np.argmax(distanser))

    print(f"Startenergi i svinghjulet ved {a.rpm:.0f} o/min: {resultater[0]['E0']:.2f} J")
    print(f"{'i':>5} {'distanse [m]':>13} {'tid [s]':>8} {'v_maks [m/s]':>13}")
    for i, r in zip(utvekslinger, resultater):
        merke = "  <- beste" if r is resultater[beste] else ""
        print(f"{i:5.1f} {r['distanse']:13.2f} {r['tid']:8.1f} {r['v_maks']:13.2f}{merke}")

    r = resultater[beste]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(utvekslinger, distanser, "o-", ms=3)
    ax1.set(xlabel="Utveksling i [-]", ylabel="Kjøredistanse [m]", title="Distanse som funksjon av utveksling")
    ax1.grid(alpha=0.3)
    ax2.plot(r["t"], r["v"])
    ax2.set(xlabel="Tid [s]", ylabel="Fart [m/s]", title=f"Fartsforløp, i = {utvekslinger[beste]:.0f}")
    ax2.grid(alpha=0.3)
    fig.savefig(a.figur, dpi=150, bbox_inches="tight")
    print(f"Figur lagret som {a.figur}")


if __name__ == "__main__":
    main()
