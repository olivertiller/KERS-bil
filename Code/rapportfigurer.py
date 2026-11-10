"""Figurer til rapporten fra måledataene for utsnurringen (ML- og TR-serien) og simuleringen.

    python3 Code/rapportfigurer.py

Skriver PDF-er til Rapport/Images/Maling/:
  ml_flanker.pdf          tid mellom flankene i ML-1 (sensoren for nær) og ML-2
  utsnurring_referanse.pdf  referanseserien (tannhjulet koblet) med tilpasset tapsmodell
  tapsmoment.pdf          tapsmoment mot turtall, tannhjulet koblet og frakoblet
  lagervarianter.pdf      tapsmoment mot turtall for fire lagervarianter, tannhjulet frakoblet
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
TORR1 = DATA / "diagnose_20261107_001343.csv"  # 626 2RS uten tetninger og uten olje, første serie (innkjøring)
TORR2 = DATA / "diagnose_20261107_002109.csv"  # samme lagre, andre serie rett etter
# diagnose_20261107_000144.csv er en prøvekjøring med samme lagre og lave opptrekk, og brukes ikke.
UTELAT_REF, UTELAT_FRA = (1, 5), (1,)
LAGERVARIANTER = (("626 2RS med fett", FRA, UTELAT_FRA), ("626 2RS uten tetninger, olje", BAMBU, ()),
                  ("626 2RS uten tetninger, tørt", TORR2, ()),
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
    ax[0].set(xlabel="Flanke nr.", ylabel="Tid til neste flanke [ms]", title="Første oppsett: sensoren for nær, 6 flanker per omdreining")
    t2 = flanker(REF)
    biter = [b for b in np.split(t2, np.where(np.diff(t2) > 1.0)[0] + 1) if len(b) > 100]
    d2 = np.diff(biter[2]) * 1e3  # trekk 3, midt i utsnurringen
    b = len(d2) // 3
    ax[1].plot(np.arange(30), d2[b:b + 30], "o-", ms=3, color="C1")
    ax[1].set(xlabel="Flanke nr.", title="Andre oppsett: sensoren flyttet, 2 flanker per omdreining")
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


def tab_lagerplan(p_ref):
    """Simulert distanse fra 2700 o/min med i = 73/17 for lagerkonfigurasjonene i tab:lagerplan.
    Tap i lagrene på en hjulaksel regnes om til ekstra rullemotstand T/r_h. Optimering av
    utvekslingen (tidligere fig_sim / sim_maalte_tap.pdf) er flyttet til KERS-bil-dynamikkV27."""
    b0 = ks.Bil(r_h=0.034, J_hjul=2.9e-5)
    p_skf, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(SKF, ())], ks.Bil().I)
    ref = replace(b0, c0=p_ref[0], c1=p_ref[1], c2=p_ref[2])
    skf = replace(b0, c0=p_skf[0], c1=p_skf[1], c2=p_skf[2])
    i, rpm = 73 / 17, 2700
    ekstra = lambda b, T: replace(b, C_rr=b.C_rr + T / b.r_h / (b.m * ks.g))
    p_torr, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(TORR2, ())], ks.Bil().I)
    torr = replace(b0, c0=p_torr[0], c1=p_torr[1], c2=p_torr[2])
    rader = (("Opprinnelige 626 2RS, tannhjulet koblet (referanseserien)", ref),
             ("SKF på svinghjulsakselen, ingen andre lagertap", skf),
             ("  + 626 2RS-par med fett foran, 1,3 mNm", ekstra(skf, 1.3e-3)),
             ("  + ett 626 2RS foran (trehjuling), 0,65 mNm", ekstra(skf, 0.65e-3)),
             ("  + SKF-par foran, 0,14 mNm", ekstra(skf, 0.14e-3)),
             ("  + ett SKF foran (trehjuling), 0,07 mNm", ekstra(skf, 0.07e-3)),
             ("  + par tørre 626 2RS uten tetninger foran, 0,6 mNm", ekstra(skf, 0.6e-3)),
             ("Tørre 626 2RS på alle akslene (kjøretesten), 0,6 mNm per hjulaksel", ekstra(torr, 1.2e-3)),
             ("SKF på svinghjulsakselen, tørre 626 2RS på hjulakslene, 0,6 mNm per par", ekstra(skf, 1.2e-3)),
             ("Tørre 626 2RS på alle akslene (kjøretesten), 1,5 mNm per hjulaksel", ekstra(torr, 3.0e-3)))
    print(f"Lagerplan, {rpm} o/min, i = 4,29:")
    for navn, b in rader:
        print(f"  {navn:60s} {ks.simuler(b, i, rpm)['distanse']:5.1f} m")
    print("Kjøretesten over turtallsområdet for snoropptrekk (1700 og 2850 o/min):")
    for navn, b in rader[7::2]:
        print(f"  {navn:60s} " + "  ".join(f"{ks.simuler(b, i, r)['distanse']:5.1f}" for r in (1700, 2850)) + " m")


def fig_maalegeometri():
    """Skisse av hvordan lengde og sideavvik måles i kjøretesten, med forsøk 4 som eksempel
    (L = 6,14 m, y = -2,80 m, R = L^2 / (2|y|))."""
    L, y = 6.14, -2.80
    R = L**2 / (2 * abs(y))
    phi_s = math.acos(1 - abs(y) / R)
    phi = np.linspace(0, phi_s, 100)
    x_b, y_b = -R * (1 - np.cos(phi)), R * np.sin(phi)  # x mot høyre, y fremover
    xs, ys = x_b[-1], y_b[-1]
    fig, ax = plt.subplots(figsize=(2.6, 3.2))
    ax.plot([0, 0], [-0.3, ys + 0.6], color="0.4", ls="--", lw=1)
    ax.text(0.12, ys + 0.45, "startlinje", fontsize=7, color="0.3")
    ax.plot(x_b, y_b, color="C0", lw=2)
    phi_t = math.asin(3.6 / R)
    ax.annotate("bilens\nbane", xy=(-R * (1 - math.cos(phi_t)), 3.6), xytext=(0.15, 2.6), fontsize=7, color="C0",
                arrowprops=dict(arrowstyle="->", color="C0", lw=0.8))
    ax.plot([0, xs], [0, ys], color="C3", lw=1.2)
    ax.text(xs / 2 - 0.95, ys / 2 - 0.3, "lengde $L$", fontsize=7, color="C3", rotation=63)
    ax.annotate("", xy=(xs, ys), xytext=(0, ys), arrowprops=dict(arrowstyle="<->", color="C1", lw=1.2))
    ax.text(xs / 2, ys + 0.18, "sideavvik $y$", fontsize=7, color="C1", ha="center")
    ax.plot(0, 0, "ks", ms=5)
    ax.text(0.15, -0.15, "start", fontsize=7, va="center")
    ax.plot(xs, ys, "ko", ms=4)
    ax.text(xs - 0.15, ys - 0.35, "stopp", fontsize=7, ha="center")
    ax.set_aspect("equal")
    ax.set_xlim(xs - 0.7, 1.3)
    ax.axis("off")
    fig.tight_layout(pad=0.1)
    ut = ROT / "Rapport" / "Images" / "Kjoretest"
    ut.mkdir(parents=True, exist_ok=True)
    fig.savefig(ut / "maalegeometri.pdf")
    plt.close(fig)


KJORETEST1 = ((3.29, -0.50), (6.12, -4.15), (6.76, -1.75), (6.14, -2.80), (5.09, 0.0))  # 7. nov, (L, y) [m]
KJORETEST2 = ((8.47, -2.47), (8.09, -0.83), (10.93, -3.91), (8.34, -1.56), (11.86, -3.24),
              (7.08, -0.19), (9.25, -2.78), (7.63, -1.12), (10.54, -3.65), (5.89, -2.04))  # 9. nov


def fig_kjoretester():
    """Banene i de to kjøretestene, tegnet som sirkelbuer fra start til stopp med R = L^2 / (2|y|).
    Det skyggelagte feltet er der |y|/L < 10 % (R5). Data fra Bygg/Testing/Testing.xlsx."""
    fig, ax = plt.subplots(figsize=(3.2, 4.6))
    vinkel = math.asin(0.10)
    yy = np.array([0, 12.5])
    ax.fill_betweenx(yy, -yy * math.tan(vinkel), yy * math.tan(vinkel), color="0.85", lw=0, zorder=0,
                     label="R5: $|y|/L$ < 10 %")
    for data, farge, navn in ((KJORETEST1, "C1", "Første kjøretest"), (KJORETEST2, "C0", "Andre kjøretest")):
        for k, (L, y) in enumerate(data):
            if y == 0:
                xb, yb = np.array([0, 0]), np.array([0, L])
            else:
                R = L**2 / (2 * abs(y))
                phi = np.linspace(0, math.acos(1 - abs(y) / R), 100)
                xb, yb = -R * (1 - np.cos(phi)), R * np.sin(phi)  # x mot høyre, y fremover
            ax.plot(xb, yb, color=farge, lw=1.5, label=navn if k == 0 else None)
            ax.plot(xb[-1], yb[-1], "o", color=farge, ms=4, mec="white", mew=0.8)
    ax.plot(0, 0, "ks", ms=5)
    ax.set_aspect("equal")
    ax.set_xlim(-4.6, 1.6)
    ax.set_ylim(-0.3, 12.5)
    ax.set_xlabel("Sideavvik $y$ [m]")
    ax.set_ylabel("Fremover [m]")
    ax.legend(loc="lower left", fontsize=7, framealpha=1)
    fig.tight_layout(pad=0.2)
    ut = ROT / "Rapport" / "Images" / "Kjoretest"
    ut.mkdir(parents=True, exist_ok=True)
    fig.savefig(ut / "kjoretester.pdf")
    plt.close(fig)


def tab_dv_versjoner():
    """Simulert distanse for drivverkversjonene i tab:dv-samlet, med dagens modell: kjøretestbilen (tørre
    626 2RS, 0,6 mNm per hjulaksel), 0,27 kg, 2700 o/min. Hjulenes treghetsmoment skaleres med r^3."""
    p, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(TORR2, ())], ks.Bil().I)
    for v, z, r_h in ((1, 68, 0.032), (2, 68, 0.034), (4, 73, 0.034)):
        b = ks.Bil(r_h=r_h, J_hjul=2.9e-5 * (r_h / 0.034) ** 3, c0=p[0], c1=p[1], c2=p[2])
        b = replace(b, C_rr=b.C_rr + 1.2e-3 / r_h / (b.m * ks.g))
        print(f"Drivverk versjon {v}: i = {z}/17, r_h = {r_h * 1e3:.0f} mm: {ks.simuler(b, z / 17, 2700)['distanse']:.2f} m")


def tab_tapsledd():
    """Energiregnskap for simulert kjøring med kjøretestkonfigurasjonen (tab:tapsledd), og følsomhet for mu."""
    b0 = ks.Bil(r_h=0.034, J_hjul=2.9e-5)
    p, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(TORR2, ())], ks.Bil().I)
    torr = replace(b0, c0=p[0], c1=p[1], c2=p[2])
    for mu in (0.6, 0.3, 1.0):
        bil = replace(torr, mu=mu, C_rr=torr.C_rr + 1.2e-3 / torr.r_h / (torr.m * ks.g))
        r = ks.simuler(bil, 73 / 17, 2700)
        E, E0 = r["E"], r["E0"]
        hjul = 1 - b0.C_rr / bil.C_rr  # andel av rullemotstanden som er hjulaksellagrene
        print(f"Tapsledd, mu = {mu}: E0 = {E0:.2f} J, {r['distanse']:.1f} m")
        for navn, e in (("Slur ved nedsetting", E["slur"] + E["las"]), ("Lagre svinghjulsaksel", E["lager"]),
                        ("Lagre hjulaksler", E["rull"] * hjul), ("Tannhjulsinngrep", E["tannhjul"]),
                        ("Rullemotstand dekk", E["rull"] * (1 - hjul))):
            print(f"  {navn:24s} {e:5.2f} J  {e / E0 * 100:4.0f} %")


def fig_sim_kjoring():
    """Simulert kjøring for kjøretestkonfigurasjonen (tørre 626 2RS, 0,6 mNm per hjulaksel), 2700 o/min.
    Viser svinghjulets turtall og bilens fart mot tid, med slurfasen ved nedsetting markert."""
    b0 = ks.Bil(r_h=0.034, J_hjul=2.9e-5)
    p, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(TORR2, ())], ks.Bil().I)
    torr = replace(b0, c0=p[0], c1=p[1], c2=p[2])
    bil = replace(torr, C_rr=torr.C_rr + 1.2e-3 / torr.r_h / (torr.m * ks.g))
    i = 73 / 17
    r = ks.simuler(bil, i, 2700)
    k = bil.r_h / i
    slur = np.abs(k * r["w"] - r["v"]) > 1e-3
    t_slur = r["t"][slur].max() if slur.any() else 0.0
    fig, ax = plt.subplots(figsize=(6.3, 2.4))
    ax.plot(r["t"], r["w"] * 30 / math.pi, color="C0", lw=1.4)
    ax.set(xlabel="Tid [s]", ylabel="Svinghjul [o/min]", xlim=(0, r["tid"]), ylim=(0, 2900))
    ax.tick_params(axis="y", colors="C0")
    ax.yaxis.label.set_color("C0")
    ax2 = ax.twinx()
    ax2.plot(r["t"], r["v"], color="C1", lw=1.4)
    ax2.set_ylabel("Bilens fart [m/s]", color="C1")
    ax2.tick_params(axis="y", colors="C1")
    ax2.set_ylim(0, r["v"].max() * 1.1)
    ax.axvspan(0, t_slur, color="0.85", lw=0)
    ax.annotate(f"slur, {t_slur:.2f} s", xy=(t_slur, 2600), xytext=(t_slur + 0.4, 2600), fontsize=8,
                va="center", arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.text(0.98, 0.9, f"{r['distanse']:.1f} m på {r['tid']:.1f} s", transform=ax.transAxes, ha="right", fontsize=8)
    fig.tight_layout()
    ut = ROT / "Rapport" / "Images" / "Simulering"
    ut.mkdir(parents=True, exist_ok=True)
    fig.savefig(ut / "sim_kjoring.pdf")
    plt.close(fig)
    print(f"Simulert kjøring: {r['distanse']:.2f} m, {r['tid']:.1f} s, v_maks {r['v_maks']:.2f} m/s, slur {t_slur:.3f} s")


def fig_lager():
    """Tapsmoment mot turtall for de fire lagervariantene på svinghjulsakselen, logaritmisk akse."""
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
    p, _ = ks.tilpass_utsnurring([sn for _, sn in trekk(TORR2, ())], ks.Bil().I)
    print(f"Tørt: c0 = {p[0] * 1e3:.3f} mNm, c1 = {p[1]:.1e}, c2 = {p[2]:.1e} N m s^2")
    print("Innkjøring, første serie med tørre lagre (trekk, topp, til stopp, 1000->300, moment):")
    for r in tid_1000_300(TORR1, ()):
        print(f"  {r[0]}  {r[1]:5.0f}  {r[2]:6.2f}  {r[3]:5.2f}  {r[4]:4.2f}")


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
    tab_lagerplan(p)
    fig_sim_kjoring()
    fig_maalegeometri()
    fig_kjoretester()
    tab_dv_versjoner()
    tab_tapsledd()
    fig_lager()
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
