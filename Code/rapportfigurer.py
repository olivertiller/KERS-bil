"""Figurer til rapporten fra måledataene for utsnurringen (ML- og TR-serien) og simuleringen.

    python3 Code/rapportfigurer.py

Skriver PDF-er til Rapport/Images/Maling/:
  ml_flanker.pdf          tid mellom flankene i ML-1 (sensoren for nær) og ML-2
  utsnurring_referanse.pdf  referanseserien (tannhjulet koblet) med tilpasset tapsmodell
  tapsmoment.pdf          tapsmoment mot turtall, tannhjulet koblet og frakoblet
  sim_maalte_tap.pdf      simulert distanse mot utveksling med antatte og målte tap
  lagervarianter.pdf      tapsmoment mot turtall for tre lagervarianter, tannhjulet frakoblet
Tallene som skrives ut, er de som står i rapporten.
"""

import math
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kers_sim as ks

ROT = Path(__file__).resolve().parent.parent
DATA = ROT / "Code" / "utsnurring_diagnose" / "data"
UT = ROT / "Rapport" / "Images" / "Maling"
ML1 = DATA / "diagnose_20261105_171409.csv"  # drill, trinn 1, sensoren for nær
REF = DATA / "diagnose_20261106_022052.csv"  # tannhjulet koblet, fett i lagrene
FRA = DATA / "diagnose_20261106_022930_frakoblet.csv"  # tannhjulet frakoblet
BAMBU = DATA / "diagnose_20261106_161745.csv"  # samme 626 2RS-lagre, tetningene ut, Bambu Lab-olje
SKF = DATA / "diagnose_20261106_163914.csv"  # SKF Explorer 626, åpne lagre
UTELAT_REF, UTELAT_FRA = (1, 5), (1,)
LAGERVARIANTER = (("626 2RS med fett", FRA, UTELAT_FRA), ("626 2RS uten tetninger, olje", BAMBU, ()),
                  ("SKF Explorer 626, åpent", SKF, ()))
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})


def flanker(fil):
    return np.loadtxt(fil, delimiter=",", skiprows=1, usecols=0, ndmin=1)


def trekk(fil, utelat, kutt=ks.KUTT_ETTER_TOPP):
    serier = ks.utsnurring_fra_flanker(flanker(fil), 2, kutt=kutt)
    return [(n, s) for n, s in enumerate(serier, 1) if n not in utelat]


def moment_mot_turtall(t, rpm, vindu=8, steg=4):
    w = rpm * math.pi / 30
    W, M = [], []
    for a in range(0, len(t) - vindu, steg):
        sl = slice(a, a + vindu)
        W.append(rpm[sl].mean())
        M.append(-np.polyfit(t[sl], w[sl], 1)[0] * ks.Bil().I * 1e3)
    return np.array(W), np.array(M)


def fig_flanker():
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.2), sharey=False)
    t1 = flanker(ML1)
    d1 = np.diff(t1) * 1e3
    a = len(d1) // 2
    ax[0].plot(np.arange(30), d1[a:a + 30], "o-", ms=3)
    ax[0].set(xlabel="Flanke nr.", ylabel="Tid til neste flanke [ms]", title="ML-1: sensoren for nær, 6 flanker per omdreining")
    t2 = flanker(REF)
    biter = [b for b in np.split(t2, np.where(np.diff(t2) > 1.0)[0] + 1) if len(b) > 100]
    d2 = np.diff(biter[2]) * 1e3  # trekk 3, midt i utsnurringen
    b = len(d2) // 3
    ax[1].plot(np.arange(30), d2[b:b + 30], "o-", ms=3, color="C1")
    ax[1].set(xlabel="Flanke nr.", title="ML-2: sensoren flyttet, 2 flanker per omdreining")
    ax[1].set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(UT / "ml_flanker.pdf")
    plt.close(fig)
    periode = d1[a:a + 60]
    print(f"ML-1: tider mellom flankene (ms), utdrag: {np.round(periode[:6], 1)}")


def fig_referanse():
    serier = trekk(REF, UTELAT_REF)
    p, w_fit = ks.tilpass_utsnurring([s for _, s in serier], ks.Bil().I)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for (n, (t, rpm)), w in zip(serier, w_fit):
        l, = ax.plot(t - t[0], rpm, ".", ms=2, label=f"trekk {n}")
        ax.plot(t - t[0], w * 30 / math.pi, color=l.get_color(), lw=0.8)
    ax.set(xlabel="Tid fra start av analysen [s]", ylabel="Turtall [o/min]",
           title="Utsnurring med tannhjulet koblet (linjer: tilpasset modell)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(UT / "utsnurring_referanse.pdf")
    plt.close(fig)
    return p


def tid_1000_300(fil, utelat):
    ut = []
    for n, (t, rpm) in trekk(fil, utelat):
        dt = ks.tid_mellom(t, rpm, 1000, 300)
        moment = ks.Bil().I * 700 * math.pi / 30 / dt * 1e3
        topp = ks.utsnurring_fra_flanker(flanker(fil), 2, kutt=0.0)[n - 1][1][0]
        ut.append((n, topp, t[-1] - t[0] + ks.KUTT_ETTER_TOPP, dt, moment))
    return ut


def fig_moment():
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    for a, fil, utelat, tittel in ((ax[0], REF, UTELAT_REF, "Tannhjulet koblet"),
                                   (ax[1], FRA, UTELAT_FRA, "Tannhjulet frakoblet")):
        for n, (t, rpm) in trekk(fil, utelat):
            W, M = moment_mot_turtall(t - t[0], rpm)
            a.plot(W, M, ".-", ms=2, lw=0.5, label=f"trekk {n}")
        a.set(xlabel="Turtall [o/min]", title=tittel, ylim=(0, 7))
        a.legend(fontsize=7)
    ax[0].set_ylabel("Tapsmoment $I\\,\\mathrm{d}\\omega/\\mathrm{d}t$ [mNm]")
    fig.tight_layout()
    fig.savefig(UT / "tapsmoment.pdf")
    plt.close(fig)


def fig_sim(p):
    b0 = ks.Bil(r_h=0.034, J_hjul=2.9e-5)
    bm = replace(b0, c0=p[0], c1=p[1], c2=p[2])
    ii = np.arange(1.5, 10.01, 0.25)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for a, b, tittel in ((ax[0], b0, "Antatte lagertap"), (ax[1], bm, "Målte tap (referanseserien)")):
        for rpm in (4000, 2700):
            d = [ks.simuler(b, x, rpm)["distanse"] for x in ii]
            k = int(np.argmax(d))
            r = ks.simuler(b, 73 / 17, rpm)["distanse"]
            a.plot(ii, d, label=f"{rpm} o/min")
            print(f"{tittel}, {rpm} o/min: i = 4,29 gir {r:.1f} m, beste i = {ii[k]:.2f} gir {d[k]:.1f} m")
        a.axvline(73 / 17, color="k", ls=":", lw=0.8)
        a.set(xlabel="Utveksling $i$", title=tittel)
        a.legend(fontsize=7)
    ax[0].set_ylabel("Simulert kjøredistanse [m]")
    fig.tight_layout()
    fig.savefig(UT / "sim_maalte_tap.pdf")
    plt.close(fig)


def fig_lager():
    """Tapsmoment mot turtall for de tre lagervariantene på svinghjulsakselen, logaritmisk akse."""
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for k, (navn, fil, utelat) in enumerate(LAGERVARIANTER):
        for j, (n, (t, rpm)) in enumerate(trekk(fil, utelat)):
            W, M = moment_mot_turtall(t - t[0], rpm, vindu=16, steg=8)
            ax.plot(W, M, ".", ms=1.5, color=f"C{k}", alpha=0.5, label=navn if j == 0 else None)
    ax.set_yscale("log")
    ax.set(xlabel="Turtall [o/min]", ylabel="Tapsmoment $I\\,\\mathrm{d}\\omega/\\mathrm{d}t$ [mNm]",
           title="Lagrene på svinghjulsakselen, tannhjulet frakoblet", ylim=(0.05, 10))
    ax.legend(fontsize=7, markerscale=6, loc="center right")
    fig.tight_layout()
    fig.savefig(UT / "lagervarianter.pdf")
    plt.close(fig)
    for navn, fil, utelat in LAGERVARIANTER[1:]:
        rader = tid_1000_300(fil, utelat)
        dt = np.array([r[3] for r in rader])
        print(f"{navn}: trekk, topp [o/min], til stopp [s], 1000->300 [s], moment [mNm]")
        for r in rader:
            print(f"  {r[0]}  {r[1]:5.0f}  {r[2]:6.2f}  {r[3]:5.2f}  {r[4]:4.2f}")
        topp = np.array([r[1] for r in rader])
        print(f"  snitt 1000->300 {dt.mean():.2f} ± {dt.std(ddof=1):.2f} s ({dt.std(ddof=1) / dt.mean() * 100:.1f} %), "
              f"moment {ks.Bil().I * 700 * math.pi / 30 / dt.mean() * 1e3:.2f} mNm, topp {topp.mean():.0f} ± {topp.std(ddof=1):.0f} o/min")
    p, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(SKF, ())], ks.Bil().I)
    print(f"SKF: c0 = {p[0] * 1e3:.3f} mNm, c1 = {p[1]:.1e}, c2 = {p[2]:.1e} N m s^2")
    b = ks.Bil(r_h=0.034, J_hjul=2.9e-5, c0=p[0], c1=p[1], c2=p[2])
    ii = np.arange(1.5, 10.01, 0.25)
    for rpm in (2700, 4000):
        d = [ks.simuler(b, x, rpm)["distanse"] for x in ii]
        k = int(np.argmax(d))
        print(f"SKF, {rpm} o/min: i = 4,29 gir {ks.simuler(b, 73 / 17, rpm)['distanse']:.1f} m, "
              f"beste i = {ii[k]:.2f} gir {d[k]:.1f} m")


def main():
    UT.mkdir(parents=True, exist_ok=True)
    fig_flanker()
    p = fig_referanse()
    print(f"Tapsmodell: c0 = {p[0] * 1e3:.2f} mNm, c1 = {p[1]:.1e}, c2 = {p[2]:.2e} N m s^2")
    for navn, fil, ut in (("Koblet", REF, UTELAT_REF), ("Frakoblet", FRA, UTELAT_FRA)):
        rader = tid_1000_300(fil, ut)
        dt = np.array([r[3] for r in rader])
        print(f"{navn}: trekk, topp [o/min], til stopp [s], 1000->300 [s], moment [mNm]")
        for r in rader:
            print(f"  {r[0]}  {r[1]:5.0f}  {r[2]:5.2f}  {r[3]:5.2f}  {r[4]:4.2f}")
        topp = np.array([r[1] for r in rader])
        print(f"  snitt 1000->300 {dt.mean():.2f} ± {dt.std(ddof=1):.2f} s ({dt.std(ddof=1) / dt.mean() * 100:.1f} %), "
              f"topp {topp.mean():.0f} ± {topp.std(ddof=1):.0f} o/min")
    fig_moment()
    fig_sim(p)
    fig_lager()
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
