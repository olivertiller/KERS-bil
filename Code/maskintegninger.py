"""Arbeids- og sammenstillingstegninger for KERS-bilen, A4 i målestokk (førstevinkelprojeksjon).

    python3 Code/maskintegninger.py            # alle tegningene
    python3 Code/maskintegninger.py TR-0 FE-1  # bare noen

Geometrien hentes fra de samme modellene som 3MF/STL-filene som printes (snortrommel.py, pinjong.py,
hjul.py, avstandshylse.py, rammen og lagerlokket fra Fusion). Riss og snitt beregnes fra modellen
med tegneverktoy.py, så kantene bak et snittplan og skjulte kanter blir riktige.
Skriver Rapport/Images/Tegninger/KERS-<nr>.pdf. Arkene er A4 i full størrelse og må settes inn i
rapporten i 100 % (\\includepdf), ellers stemmer ikke målestokken.
"""

import math
import sys
from pathlib import Path

import manifold3d as mf
import numpy as np

import avstandshylse as AH
import hjul as H
import pinjong as PI
import snortrommel as T
import tegneverktoy as tv

sys.path.insert(0, str(Path(__file__).resolve().parent / "Fusion" / "KERS_GD_oppsett"))
from tegneverktoy import Ark, Riss, tall

ROT = Path(__file__).resolve().parent.parent
UT = ROT / "Rapport" / "Images" / "Tegninger"
X, Y, Z = np.eye(3)


def masse(m, rho=T.RHO_PLA):
    return f"{tall(round(m.volume() * rho, 1))} g"


# ==========================================================================================
# Snortromler TR-0 og TR-1
# ==========================================================================================

def tegn_trommel(variant):
    tr0 = variant == "TR-0"
    p, d_k = (T.P_TR0, 8.0) if tr0 else (T.P, 12.0)
    m = T.trommel(p, d_k)
    rk, rf, (z1, z2, z3, z4, z5) = T.lengder(p, d_k)
    rb, rn, ru = p["boring"] / 2, p["d_nav"] / 2, p["d_utsparing"] / 2
    ri = p.get("d_indre", 2 * rf) / 2
    zm = (z4 + z5) / 2
    zs = z1 + p["krok_klaring"]
    ze = zs + p["krok_bredde"]
    e, xs = p["e_skrue"], p["x_sete"]
    k = 2.0
    ark = Ark()

    # --- Snitt A-A (hovedriss): planet y = 0, aksen vannrett, navet til høyre, kroken ned
    A = Riss(ark, m, Z, X, (128, 148), k, snitt=(0, 0, 0))
    a = lambda z, x: A.p(x, 0, z)
    A.senterlinje((0, 0, -3), (0, 0, z5 + 3))
    # aksiale mål fra flaten mot svinghjulet (z = 0), parallellmål uten kjedemål
    yb = a(0, rn)[1] + 3
    if tr0:
        trinn = [(p["rand_indre"], ri), (z1, rk), (z2, rk), (z5, rn)]
    else:
        trinn = [(z1, rf), (z2, rk), (z4, rf), (z5, rn)]
    for i, (zz, rr) in enumerate(trinn):
        ark.maal_h(a(0, ri), a(zz, rr), yb + 5 * i, tall(zz))
    ark.overskrift(a(z5 / 2, 0)[0], yb + 5 * len(trinn) + 1, "A–A")
    # diametre
    ark.maal_v(a(0, -ri), a(0, ri), a(-9, 0)[0], "Ø" + tall(2 * ri))
    ark.maal_v(a(0, -ru), a(0, ru), a(-4.5, 0)[0], "Ø" + tall(2 * ru))
    ark.maal_v(a(z5, -rb), a(z5, rb), a(z5 + 3.5, 0)[0], "Ø" + tall(2 * rb))
    ark.maal_v(a(z5, -rn), a(z5, rn), a(z5 + 8, 0)[0], "Ø" + tall(2 * rn))
    # utsparingen: rett del og (TR-0) flatt tak
    zr = p["rett_utsparing"]
    yu = a(0, -ri)[1] - 4
    ark.maal_h(a(0, -ru), a(zr, -ru), yu, tall(zr), ute=-1)
    if tr0:
        zt = zr + ru - rk
        ark.maal_h(a(0, -ri), a(zt, -rk), yu - 5, tall(zt), ute=1)
    # 45°-vinkler, målt på luftsiden mot en loddrett hjelpelinje
    c = a(z2, rk)  # kjegla fra kjernen ut mot flensen/navet
    ark.vinkelmaal(c, 8, 45, 90, "45°", hjelpelinjer=[(c + [0, 1], c + [0, 10])])
    if tr0:
        c2 = a(z1, rk)  # indre flanke
        ark.vinkelmaal(c2, 8, 90, 135, "45°", hjelpelinjer=[(c2 + [0, 1], c2 + [0, 10])])
    cu = a(zr, ru)
    ark.vinkelmaal(cu, 9, 270, 315, "45°", hjelpelinjer=[(cu - [0, 1], cu - [0, 11])])
    # kroken: aksial plassering og bredde (målt i snittet, undersiden)
    yk = a(0, -rn)[1] - (14 if tr0 else 9)
    ark.maal_h(a(z1, -rk), a(zs, -rk), yk, tall(zs - z1), ute=-1)
    ark.maal_h(a(zs, -rk), a(ze, -rk), yk, tall(ze - zs))

    # --- Riss fra høyre (navenden), til venstre for hovedrisset. Snittlinje A-A.
    E_ = Riss(ark, m, -Y, X, (58, 148), k)
    e0 = E_.p(0, 0, 0)
    ark.kryss(e0, k * rn)
    ark.snittlinje(e0 + [0, k * rn + 2], [1, 0], "A", [0, 1])
    ark.snittlinje(e0 - [0, k * rn + 2], [1, 0], "A", [0, -1])
    ark.henvisning(E_.p(0.25, -(rb + 4.5), z5), e0 + [k * rn + 2, k * rn + 4], f"2 × slisse {tall(p['slisse'])}")
    ark.maal_dia(e0, k * rb, -40, "Ø" + tall(2 * rb), ute=13)

    # --- Riss ovenfra (sett fra +x), under hovedrisset: hodesete og mutterlomme
    O = Riss(ark, m, Z, -Y, (128, 76), k)
    o = lambda z, y: O.p(0, y, z)
    O.senterlinje((0, 0, -3), (0, 0, z5 + 3))
    for ys in (e, -e):
        ark.senterlinje(o(zm, ys) - [0, 5], o(zm, ys) + [0, 5])
    ark.maal_h(o(zm, -e), o(z5, -rn), o(0, -rn)[1] + 13, tall(z5 - zm))
    ark.maal_v(o(zm, e), o(zm, -e), o(z5 + 6, 0)[0], tall(2 * e))
    ark.henvisning(O.p(0, e + 2.2, zm + 1.3), o(z5 + 13, e + 5), f"Mutterlomme {tall(p['mutter'])}")
    ark.henvisning(O.p(0, -e - 1.6, zm + 2.6), o(z5 + 13, -e - 5), f"Hodesete Ø{tall(p['d_hode'])}")
    # snittlinjer B-B og C-C
    zb = 7.5 if not tr0 else ze - 0.6
    for zz, bokstav in ((zb, "B"), (zm, "C")):
        ark.snittplan(o(zz, rn + 1), o(zz, -rn - 1), [-1, 0], bokstav)

    # --- Snitt C-C gjennom skrueaksene
    C = Riss(ark, m, -Y, X, (58, 76), k, snitt=(0, 0, zm))
    c0 = C.p(0, 0, 0)
    ark.kryss(c0, k * rn)
    for ys in (e, -e):
        ark.senterlinje(C.p(-rn - 1.5, ys, zm), C.p(rn + 1.5, ys, zm))
    ark.overskrift(*(c0 + [0, k * rn + 9]), "C–C")
    ark.maal_h(C.p(rn, e, zm), C.p(rn, -e, zm), c0[1] + k * rn + 3, tall(2 * e))
    ark.maal_v(C.p(0, e, zm), C.p(-xs, e, zm), c0[0] - k * rn - 4, tall(xs))
    ark.maal_v(C.p(0, -e, zm), C.p(xs, -e, zm), c0[0] + k * rn + 4, tall(xs))
    ark.maal_h(C.p(-rn, e - p["d_hode"] / 2, zm), C.p(-rn, e + p["d_hode"] / 2, zm), c0[1] - k * rn - 5,
               "Ø" + tall(p["d_hode"]))
    ark.henvisning(C.p(-1.0, -e + p["m3_hull"] / 2, zm), c0 + [k * rn + 3, -k * rn - 5], "Ø" + tall(p["m3_hull"]))

    # --- Snitt B-B gjennom kroken, og detalj E (5:1)
    rmax = max(ri, rf)
    Bb = Riss(ark, m, -Y, X, (222, 150), k, snitt=(0, 0, zb))
    b0 = Bb.p(0, 0, 0)
    ark.kryss(b0, k * rmax)
    ark.overskrift(*(b0 + [0, k * rmax + 11]), "B–B")
    ark.maal_v(Bb.p(rk, 0, zb), Bb.p(-rk, 0, zb), b0[0] + k * rmax + 5, "Ø" + tall(2 * rk))
    # opptrekksretning: mot klokka sett fra navenden
    th = np.linspace(math.radians(25), math.radians(75), 30)
    rr = k * rmax + 4
    bue = np.c_[b0[0] + rr * np.cos(th), b0[1] + rr * np.sin(th)]
    ark.linje(bue, lw=tv.TYNN)
    ark.pil(bue[-1], bue[-1] - bue[-3])
    ark.tekst(b0[0] + rr * 0.93 + 3, b0[1] + rr * 0.36, "Opptrekk", ha="left")
    # detalj E rundt kroken (på -x, ned i risset)
    hk = p["tapp_hoyde"] if tr0 else p["krok_hoyde"]
    cE_mod = (-(rk + hk / 2), 0.6 if not tr0 else 0.3, zb)
    ark.detaljmerke(Bb.p(cE_mod), k * 3.4, "E", vinkel=210)
    kd = 5.0
    cE = np.array([248.0, 86.0])
    D = Riss(ark, m, -Y, X, cE - kd * (tv.basis(-Y, X)[:2] @ np.asarray(cE_mod)), kd, snitt=(0, 0, zb), tegn=False)
    D.tegn(klipp=(cE, kd * 4.2))
    ark.overskrift(cE[0], cE[1] + kd * 4.2 + 3, "E (5:1)")
    d = lambda t, r: D.p(-(rk + r), t, zb)  # (t langs omkretsen, r over kjerneflaten) -> arket
    if tr0:
        s, h_, ang = p["krok_stilk"] / 2, p["tapp_hoyde"], p["tapp_helning"]
        kk = math.tan(math.radians(ang))
        ark.maal_v(d(-s - h_ * kk, h_), d(-s, 0), d(-s - h_ * kk - 2.2, 0)[0], tall(h_))
        ark.vinkelmaal(d(s, 0), kd * 2.6, 90, 90 + ang, tall(ang) + "°", hjelpelinjer=[(d(s, 0), d(s, 3.0))])
        ark.henvisning(d(-s - h_ * kk + 0.2, h_ - 0.1), cE + [-14, -16], f"2 × R{tall(p['tapp_runding'])}",
                       ha="right")
        ark.maal_h(d(-s, -0.0), d(s, 0), d(0, -1.2)[1], tall(2 * s))
    else:
        s, L, g, h_ = p["krok_stilk"] / 2, p["krok_finger"], p["krok_gap"], p["krok_hoyde"]
        yb2 = d(0, h_ + 1.0)[1]
        ark.maal_h(d(-s - L, h_), d(-s, g), yb2, tall(L))
        ark.maal_h(d(-s, g), d(s, h_), yb2, tall(2 * s))
        ark.maal_v(d(s, 0), d(s, h_), d(s + 1.6, 0)[0], tall(h_))
        ark.maal_v(d(-s - L, 0), d(-s - L, g), d(-s - L - 1.3, 0)[0], tall(g))
        ark.henvisning(d(-s - L * 0.5, g * 0.45), cE + [-12, -17], "Løkka", ha="right", prikk=True)

    merk = ["Printes stående med flensen mot svinghjulet ned,", "uten support.",
            "Klemnavet strammes med 2 × M3 × 12 (ISO 7045)", "og mutter M3 (ISO 4032).",
            "Løkka på snora legges rundt kroken og glir av", "av seg selv når snora er trukket helt ut."]
    if not tr0:
        merk += ["TR-2 og TR-3: som TR-1 med kjerne Ø14 og Ø18,", "flens Ø20 og Ø24."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    tittel = "Snortrommel TR-0, kjerne Ø8" if tr0 else "Snortrommel TR-1, kjerne Ø12"
    ark.tittelfelt(tittel, f"KERS-{variant}", "PLA", masse(m), "2:1 (E 5:1)")
    ark.lagre(UT / f"KERS-{variant}.pdf")


# ==========================================================================================
# Drivhjul: felg (FE-1) og dekk (DE-1)
# ==========================================================================================

def tegn_felg():
    p = H.P
    m = H.hjul(p, True)
    b = p["bredde"]
    r_sete = p["d_ytre"] / 2 - p["t_dekk"]
    rs = r_sete + p["skulder"]
    r_inn = r_sete - p["t_felg"]
    kant = (b - p["b_dekk"]) / 2
    rn, rkr, rr_ = p["d_nav"] / 2, p["d_krage"] / 2, p["d_rund"] / 2
    zt = b + p["l_krage"]
    s2 = p["firkant"] / 2
    k = 1.0
    ark = Ark()

    # Riss sett fra utsiden (betrakteren på -z), eikene vannrett og på ±60°
    F = Riss(ark, m, X, -Y, (78, 128), k)
    f0 = F.p(0, 0, 0)
    ark.kryss(f0, k * rs)
    for v in (60, 120):
        u = np.array([math.cos(math.radians(v)), math.sin(math.radians(v))])
        ark.senterlinje(f0 - u * (k * rs + 3), f0 + u * (k * rs + 3))
    ark.snittlinje(f0 + [0, k * rs + 3], [1, 0], "A", [0, 1])
    ark.snittlinje(f0 - [0, k * rs + 3], [1, 0], "A", [0, -1])
    u = np.array([math.cos(math.radians(60)), math.sin(math.radians(60))])
    nn = np.array([-u[1], u[0]])
    pa, pb = f0 + u * k * 20 - nn * k * p["b_eike"] / 2, f0 + u * k * 20 + nn * k * p["b_eike"] / 2
    ark.henvisning(f0 + u * k * 24, f0 + [-22, 46], f"{p['n_eiker']} eiker × 60°", ha="right", prikk=True)
    ark.henvisning(f0 + k * np.array([24, 2.0]), f0 + [24, 44], f"Eikebredde {tall(p['b_eike'])}", prikk=False)
    ark.maal_dia(f0, k * r_inn, -30, "Ø" + tall(2 * r_inn))
    ark.detaljmerke(f0, k * rn + 2.5, "C", vinkel=-50)

    # Snitt A-A (planet x = 0, mellom eikene), sett mot +x. Utsiden (z = 0) til høyre.
    A = Riss(ark, m, -Z, -Y, (162, 128), k, snitt=(0, 0, 0))
    a = lambda z, r: A.p(0, -r, z)  # r > 0 oppover
    A.senterlinje((0, 0, zt + 3), (0, 0, -3))
    ark.overskrift(a(b / 2, 0)[0], a(0, rs)[1] + 14, "A–A")
    ark.maal_h(a(0, rs), a(b, rs), a(0, rs)[1] + 5, tall(b))
    ark.maal_v(a(0, -rs), a(0, rs), a(-6, 0)[0], "Ø" + tall(2 * rs))
    ark.maal_h(a(0, -r_inn + 4), a(p["t_skive"], -r_inn + 4), a(0, -rs)[1] - 5, tall(p["t_skive"]), ute=-1)
    ark.detaljmerke(a(b / 2, r_sete + 0.2), 6.0, "D", vinkel=160)
    ark.detaljmerke(a(zt / 2 + 0.4, (rn + rr_) / 2 - 0.4), 6.5, "F", vinkel=200)

    # Detalj D (5:1): felgsetet i snitt
    kd = 5.0
    cD_m = (0, -(r_sete - 0.4), b / 2)
    cD = np.array([212.0, 163.0])
    D = Riss(ark, m, -Z, -Y, cD - kd * (tv.basis(-Z, -Y)[:2] @ np.asarray(cD_m)), kd, snitt=(0, 0, 0), tegn=False)
    D.tegn(klipp=(cD, kd * 4.4))
    ark.overskrift(cD[0], cD[1] + kd * 4.4 + 3, "D (5:1)")
    d = lambda z, r: D.p(0, -r, z)
    ark.maal_h(d(b, rs), d(b - kant, rs), d(0, rs)[1] + 4, tall(kant))
    ark.maal_h(d(b - kant, r_sete), d(kant, r_sete), d(0, rs)[1] + 4, tall(b - 2 * kant))
    ark.maal_v(d(0, r_sete), d(0, rs), d(-1.2, 0)[0], tall(p["skulder"]))
    ark.maal_v(d(b, r_inn), d(b, r_sete), d(b + 1.0, 0)[0], tall(p["t_felg"]))
    ark.henvisning(d(b / 2, r_sete), cD + [26, -16], f"Felgsete Ø{tall(2 * r_sete)}", prikk=True)

    # Detalj F (5:1): navet med krage og rund forsenkning
    cF_m = (0, -4.0, zt / 2 + 0.4)
    cF = np.array([256.0, 100.0])
    Fd = Riss(ark, m, -Z, -Y, cF - kd * (tv.basis(-Z, -Y)[:2] @ np.asarray(cF_m)), kd, snitt=(0, 0, 0), tegn=False)
    Fd.tegn(klipp=(cF, kd * 5.4))
    ark.overskrift(cF[0], cF[1] + kd * 5.4 + 3, "F (5:1)")
    ark.henvisning(Fd.p(0, -rn, 1.5), cF + [-30, 22], f"Nav Ø{tall(2 * rn)}", ha="right")
    q = lambda z, r: Fd.p(0, -r, z)
    ark.maal_h(q(b, rkr), q(zt, rkr), q(0, rn)[1] + 3, tall(p["l_krage"]), ute=-1)
    ark.maal_h(q(zt - p["l_rund"], rr_), q(zt, rr_), q(0, rn)[1] + 9, tall(p["l_rund"]), ute=-1)
    ark.henvisning(q(zt - 0.5, rr_), cF + [-30, -18], f"Ø{tall(2 * rr_)}", ha="right")
    ark.henvisning(q(zt, rkr - 0.15), cF + [-30, -10], f"Ø{tall(2 * rkr)}", ha="right")

    # Detalj C (5:1): firkanthullet sett fra utsiden
    kc = 5.0
    cC = np.array([190.0, 76.0])
    Cd = Riss(ark, m, X, -Y, cC, kc, tegn=False)
    Cd.tegn(klipp=(cC, kc * 3.6))
    ark.overskrift(cC[0], cC[1] + kc * 3.6 + 3, "C (5:1)")
    ark.kryss(cC, kc * s2 + 1)
    ark.maal_h(Cd.p(-s2, s2, 0), Cd.p(s2, s2, 0), cC[1] - kc * s2 - 6, tall(p["firkant"]))
    ark.maal_v(Cd.p(s2, -s2, 0), Cd.p(s2, s2, 0), cC[0] + kc * s2 + 6, tall(p["firkant"]))
    c = s2 - 0.45 * 0.35
    ark.henvisning(Cd.p(c + 0.45 * 0.71, c + 0.45 * 0.71, 0), cC + [22, 18], "4 × R0,45")

    merk = ["Printes liggende med utsiden ned.",
            "Firkanthullet passer akselenden 4 × 4",
            "(valgt etter prøvebiten).",
            "Kragen ligger an mot innerringen i lageret.",
            "Dekket KERS-DE-1 strekkes over felgsetet."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Felg Ø65 for TPU-dekk", "KERS-FE-1", "PLA", masse(m), "1:1 (C, D, F 5:1)")
    ark.lagre(UT / "KERS-FE-1.pdf")


def tegn_dekk():
    p = H.P
    m = H.dekk(p)
    r_sete = p["d_ytre"] / 2 - p["t_dekk"]
    ri = r_sete - p["strekk"] / 2
    ro = ri + p["t_dekk"]
    bd, f = p["b_dekk"] - 0.2, 0.4
    k = 1.0
    ark = Ark()
    F = Riss(ark, m, X, -Y, (80, 125), k)
    f0 = F.p(0, 0, 0)
    ark.kryss(f0, k * ro)
    ark.snittlinje(f0 + [0, k * ro + 3], [1, 0], "A", [0, 1])
    ark.snittlinje(f0 - [0, k * ro + 3], [1, 0], "A", [0, -1])
    ark.maal_dia(f0, k * ri, 35, "Ø" + tall(2 * ri))
    ark.maal_dia(f0, k * ro, -35, "Ø" + tall(2 * ro), ute=10)
    A = Riss(ark, m, -Z, -Y, (160, 125), k, snitt=(0, 0, 0))
    a = lambda z, r: A.p(0, -r, z)
    A.senterlinje((0, 0, bd + 3), (0, 0, -3))
    ark.overskrift(a(bd / 2, 0)[0], a(0, ro)[1] + 12, "A–A")
    ark.maal_h(a(0, ro), a(bd, ro), a(0, ro)[1] + 5, tall(bd))
    ark.detaljmerke(a(bd / 2, ro - 0.5), 6.0, "B", vinkel=200)
    kd = 5.0
    cB_m = (0, -(ri + p["t_dekk"] / 2), bd / 2)
    cB = np.array([225.0, 140.0])
    D = Riss(ark, m, -Z, -Y, cB - kd * (tv.basis(-Z, -Y)[:2] @ np.asarray(cB_m)), kd, snitt=(0, 0, 0), tegn=False)
    D.tegn(klipp=(cB, kd * 4.4))
    ark.overskrift(cB[0], cB[1] + kd * 4.4 + 3, "B (5:1)")
    d = lambda z, r: D.p(0, -r, z)
    ark.maal_h(d(0, ri), d(bd, ri), d(0, ri)[1] - 6, tall(bd))
    ark.maal_v(d(bd, ri), d(bd, ro), d(bd + 1.4, 0)[0], tall(p["t_dekk"]))
    ark.henvisning(d(bd - f / 2, ro - f / 2), cB + [16, 14], f"2 × {tall(f)} × 45°")
    merk = ["Strekkes over felgsetet Ø65 på KERS-FE-1",
            f"(strekk {tall(p['strekk'])} mm på diameteren). Ferdig hjul Ø68.",
            "Printes liggende i TPU, uten support."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Dekk for drivhjul", "KERS-DE-1", "TPU", masse(m, H.RHO_TPU), "1:1 (B 5:1)")
    ark.lagre(UT / "KERS-DE-1.pdf")


# ==========================================================================================
# Avstandshylse AH-1
# ==========================================================================================

def tegn_hylse():
    p = AH.P
    L = 3.0
    m = AH.hylse(p, L)
    lk = p["l_krage"]
    rb, ry, rk = p["boring"] / 2, p["d_ytre"] / 2, p["d_krage"] / 2
    k = 5.0
    ark = Ark()
    # Riss sett fra kragen (betrakteren på +z), åpningen mot høyre
    E = Riss(ark, m, X, Y, (75, 128), k)
    e0 = E.p(0, 0, 0)
    ark.kryss(e0, k * ry)
    ark.maal_dia(e0, k * rb, 125, "Ø" + tall(2 * rb), ute=12)
    ark.maal_dia(e0, k * ry, 215, "Ø" + tall(2 * ry), ute=8)
    a1, a2 = p["apning"] / 2, p["apning_ytre"] / 2
    ark.maal_v(E.p(rb, -a1, 0), E.p(rb, a1, 0), e0[0] + k * ry + 6, tall(p["apning"]))
    ark.maal_v(E.p(ry, -a2, 0), E.p(ry, a2, 0), e0[0] + k * ry + 16, tall(p["apning_ytre"]))
    ark.snittlinje(e0 + [-k * ry - 3, 0], [0, 1], "A", [-1, 0])
    ark.snittlinje(e0 + [k * ry + 3, 0], [0, 1], "A", [1, 0])
    # Snitt A-A (planet y = 0, gjennom åpningen), sett mot -y, under risset
    A = Riss(ark, m, X, Z, (75, 52), k, snitt=(0, 0, 0))
    ark.senterlinje(A.p(0, 0, -2) , A.p(0, 0, L + 2))
    ark.overskrift(A.p(0, 0, 0)[0], A.p(0, 0, L)[1] + 14, "A–A")
    ark.maal_v(A.p(-ry, 0, 0), A.p(-ry, 0, L), A.p(-ry, 0, 0)[0] - 7, "L")
    ark.maal_v(A.p(rk, 0, L - lk), A.p(rk, 0, L), A.p(ry, 0, 0)[0] + 7, tall(lk))
    ark.maal_h(A.p(-rk, 0, L), A.p(rk, 0, L), A.p(0, 0, L)[1] + 6, "Ø" + tall(2 * rk))
    # tabell over lengder
    rad = [("L [mm]", "1,0 – 7,0 i trinn på 0,5"), ("Merking", "hylse_L<L>.3mf"),
           ("Tegnet med", f"L = {tall(L)}")]
    ark.tabell(175, 150, rad, (30, 60), tittel="Lengder")
    merk = ["C-formet, klikkes på Ø6-akselen uten å ta den ut av lagrene.",
            "Kragen ligger an mot innerringen i lageret, Ø12 mot navet.",
            "L velges så akselen får 0,1–0,3 mm klaring sideveis.",
            "Printes stående med Ø12-enden ned, uten support."]
    ark.merknader(160, 110, merk, h=2.3)
    ark.tittelfelt("Avstandshylse", "KERS-AH-1", "PLA", masse(m) + " (L = 3)", "5:1")
    ark.lagre(UT / "KERS-AH-1.pdf")


# ==========================================================================================
# Lagerlokk LL-1 (fra Fusion, samme mål som KERS_GD_oppsett)
# ==========================================================================================

def tegn_lokk():
    m = tv.les_stl(ROT / "Bygg" / "Lagerlokk -y.stl").translate([-70, 0, -34])  # boringens sentrum i origo
    # bilens koordinater: x fremover, y = -25,5 (innsiden) til -32,5 (utsiden), z opp
    xb, y0, y1, z0, z1 = 16.0, -32.5, -25.5, 6.0, 13.0
    rb, rl = 19.1 / 2, 17.0 / 2
    xs, ys = 12.75, -29.5
    k = 2.0
    ark = Ark()
    # Hovedriss sett fra utsiden (betrakteren på -y)
    F = Riss(ark, m, X, Z, (115, 135), k, ref=(0, y0, 9.5), skjulte=True)
    f = lambda x, z: F.p(x, y0, z)
    ark.senterlinje(f(0, -3), f(0, z1 + 3))
    ark.senterlinje(f(-rb - 3, 0), f(rb + 3, 0))
    for x in (-xs, xs):
        ark.senterlinje(f(x, z0 - 3), f(x, z1 + 3))
    ark.maal_h(f(-xb, z1), f(xb, z1), f(0, z1)[1] + 12, tall(2 * xb))
    ark.maal_h(f(-xs, z1), f(xs, z1), f(0, z1)[1] + 6, tall(2 * xs))
    ark.maal_v(f(xb, z0), f(xb, z1), f(xb, 0)[0] + 6, tall(z1 - z0))
    ark.maal_v(f(rb + 3, 0), f(xb, z1), f(xb, 0)[0] + 15, tall(z1))
    ark.henvisning(f(rb * math.cos(math.radians(120)), rb * math.sin(math.radians(120))),
                   f(-xb - 6, z0 - 6), f"Ø{tall(2 * rb)}", ha="right")
    ark.snittplan(f(0, -4), f(0, z1 + 9), [1, 0], "A")
    # Riss ovenfra, under hovedrisset
    O = Riss(ark, m, X, Y, (115, 80), k, ref=(0, ys, z1))
    o = lambda x, y: O.p(x, y, z1)
    for x in (-xs, xs):
        ark.kryss(o(x, ys), k * 1.7, ut=2)
    ark.maal_v(o(xb, y0), o(xb, y1), o(xb, 0)[0] + 6, tall(y1 - y0))
    ark.maal_v(o(xs, y0), o(xs, ys), o(xb, 0)[0] + 14, tall(ys - y0))
    ark.henvisning(o(-xs - 1.2, ys + 1.2), o(-xb - 6, y1 + 6), "2 × Ø3,4 gjennom", ha="right")
    # Snitt A-A gjennom boringen (x = 0), sett mot +x, til høyre
    A = Riss(ark, m, -Y, Z, (215, 135), k, ref=(0, -29, 9.5), snitt=(0, 0, 0))
    a = lambda y, z: A.p(0, y, z)
    ark.overskrift(a(-29, z1)[0], a(0, z1)[1] + 10, "A–A")
    ark.maal_h(a(y1, z1), a(y1 - 1, z1), a(0, z1)[1] + 4, "1")
    ark.henvisning(a(y1 - 0.5, rl), a(y1, z0) + [-10, -14], f"Ø{tall(2 * rl)}", ha="right")
    ark.henvisning(a(y0 + 2, rb * 0.99), a(y0, z0) + [8, -14], f"Ø{tall(2 * rb)}")
    ark.senterlinje(a(y0 - 2, 0), a(y1 + 2, 0))
    merk = ["2 stk. Lokket på +y-siden er det samme,", "snudd 180°.",
            "Festes med 2 × M3 × 16 (ISO 7045) rett i PLA.",
            "Boring og leppe som lagerhuset i rammen", "(KERS-RA-1). Printes med oversiden ned."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Lagerlokk", "KERS-LL-1", "PLA", masse(m), "2:1")
    ark.lagre(UT / "KERS-LL-1.pdf")


# ==========================================================================================
# Pinjong PI-1 (forenklet fremstilling av tennene etter ISO 2203)
# ==========================================================================================

def klemnav_riss(ark, m, p, z4, z5, o_ende, o_topp, o_cc, k, snitt_topp=(), etiketter_venstre=False):
    """Riss fra navenden (med snittlinje A-A), riss ovenfra (med snittlinjer) og snitt C-C for et
    klemnav på z4..z5. Modellen har aksen langs z, skruene langs x på y = ±e."""
    rb, rn, e, xs = p["boring"] / 2, p["d_nav"] / 2, p["e_skrue"], p["x_sete"]
    zm = (z4 + z5) / 2
    E_ = Riss(ark, m, -Y, X, o_ende, k)
    e0 = E_.p(0, 0, 0)
    ark.kryss(e0, k * rn)
    ark.snittlinje(e0 + [0, k * rn + 2], [1, 0], "A", [0, 1])
    ark.snittlinje(e0 - [0, k * rn + 2], [1, 0], "A", [0, -1])
    ark.henvisning(E_.p(0.25, -(rb + 4.5), z5), e0 + [k * rn + 2, k * rn + 4], f"2 × slisse {tall(p['slisse'])}")
    ark.maal_dia(e0, k * rb, -40, "Ø" + tall(2 * rb), ute=13)
    O = Riss(ark, m, Z, -Y, o_topp, k)
    o = lambda z, y: O.p(0, y, z)
    O.senterlinje((0, 0, -3), (0, 0, z5 + 3))
    for ys in (e, -e):
        ark.senterlinje(o(zm, ys) - [0, 5], o(zm, ys) + [0, 5])
    ark.maal_h(o(zm, -e), o(z5, -rn), o(0, -rn)[1] + 13, tall(z5 - zm))
    ark.maal_v(o(zm, e), o(zm, -e), o(z5 + 6, 0)[0], tall(2 * e))
    if etiketter_venstre:
        ark.henvisning(O.p(0, e + 2.2, zm - 1.3), o(-5, e + 6), f"Mutterlomme {tall(p['mutter'])}", ha="right")
        ark.henvisning(O.p(0, -e - 1.6, zm - 2.6), o(-5, -e - 6), f"Hodesete Ø{tall(p['d_hode'])}", ha="right")
    else:
        ark.henvisning(O.p(0, e + 2.2, zm + 1.3), o(z5 + 13, e + 5), f"Mutterlomme {tall(p['mutter'])}")
        ark.henvisning(O.p(0, -e - 1.6, zm + 2.6), o(z5 + 13, -e - 5), f"Hodesete Ø{tall(p['d_hode'])}")
    for zz, bokstav in ((zm, "C"),) + tuple(snitt_topp):
        ark.snittplan(o(zz, rn + 1), o(zz, -rn - 1), [-1, 0], bokstav)
    C = Riss(ark, m, -Y, X, o_cc, k, snitt=(0, 0, zm))
    c0 = C.p(0, 0, 0)
    ark.kryss(c0, k * rn)
    for ys in (e, -e):
        ark.senterlinje(C.p(-rn - 1.5, ys, zm), C.p(rn + 1.5, ys, zm))
    ark.overskrift(*(c0 + [0, k * rn + 9]), "C–C")
    ark.maal_h(C.p(rn, e, zm), C.p(rn, -e, zm), c0[1] + k * rn + 3, tall(2 * e))
    ark.maal_v(C.p(0, e, zm), C.p(-xs, e, zm), c0[0] - k * rn - 4, tall(xs))
    ark.maal_v(C.p(0, -e, zm), C.p(xs, -e, zm), c0[0] + k * rn + 4, tall(xs))
    ark.maal_h(C.p(-rn, e - p["d_hode"] / 2, zm), C.p(-rn, e + p["d_hode"] / 2, zm), c0[1] - k * rn - 5,
               "Ø" + tall(p["d_hode"]))
    ark.henvisning(C.p(-1.0, -e + p["m3_hull"] / 2, zm), c0 + [k * rn + 3, -k * rn - 5], "Ø" + tall(p["m3_hull"]))


def pinjong_forenklet(p):
    """Pinjongen med tennene erstattet av toppsylinderen (forenklet fremstilling)."""
    _, (r, rg, ra, rf) = PI.tannprofil(p)
    rb = p["boring"] / 2
    z4, z5 = p["b"], p["b"] + p["l_nav"]
    m = T.syl(ra, 0, p["b"]) + T.syl(p["d_nav"] / 2, z4, z5)
    m = m - T.syl(rb, -1, z5 + 1)
    m = m - T.syl(rb + 0.4, -0.01, 0.4, rb) - T.syl(rb, z5 - 0.4, z5 + 0.01, rb + 0.4)
    return m - T.klemnav_kutt(p, rb, z4, z5), (r, rg, ra, rf)


def tegn_pinjong():
    p = PI.P
    m, (r, rg, ra, rf) = pinjong_forenklet(p)
    ekte, _ = PI.pinjong(p)
    b = p["b"]
    z4, z5 = b, b + p["l_nav"]
    rb, rn = p["boring"] / 2, p["d_nav"] / 2
    k = 2.0
    ark = Ark()
    # Snitt A-A (hovedriss)
    A = Riss(ark, m, Z, X, (138, 148), k, snitt=(0, 0, 0))
    a = lambda z, x: A.p(x, 0, z)
    A.senterlinje((0, 0, -3), (0, 0, z5 + 3))
    for sg in (1, -1):
        ark.linje([a(-2.5, sg * r), a(b + 2.5, sg * r)], lw=tv.TYNN, ls=tv.STREKPUNKT, z=5)  # delelinje
    yb = a(0, rn)[1] + 3
    ark.maal_h(a(0, ra), a(b, ra), yb, tall(b))
    ark.maal_h(a(0, ra), a(z5, rn), yb + 5, tall(z5))
    ark.overskrift(a(z5 / 2, 0)[0], yb + 11, "A–A")
    ark.maal_v(a(0, -ra), a(0, ra), a(-5, 0)[0], "Ø" + tall(2 * ra))
    ark.maal_v(a(z5, -rb), a(z5, rb), a(z5 + 3.5, 0)[0], "Ø" + tall(2 * rb))
    ark.maal_v(a(z5, -rn), a(z5, rn), a(z5 + 8, 0)[0], "Ø" + tall(2 * rn))
    ark.henvisning(a(0.2, -rb - 0.2), a(-6, -rn - 3), "2 × 0,4 × 45°", ha="right")
    # klemnavet: riss fra høyre, ovenfra og C-C
    klemnav_riss(ark, m, p, z4, z5, (58, 148), (138, 76), (52, 76), k, etiketter_venstre=True)
    # Riss fra venstre (tannsiden), til høyre for hovedrisset: forenklet fremstilling
    V = Riss(ark, m, Y, X, (215, 148), k)
    v0 = V.p(0, 0, 0)
    ark.kryss(v0, k * ra)
    ark.sirkel(v0, k * r, lw=tv.TYNN, ls=tv.STREKPUNKT)
    ark.sirkel(v0, k * rf, lw=tv.TYNN)
    # tannhjulsdata
    rad = [("Antall tenner", "z", f"{p['z']}"), ("Modul", "m", tall(p["modul"])),
           ("Trykkvinkel", "α", f"{tall(p['alfa'])}°"), ("Profilforskyvning", "x", "0"),
           ("Delesirkel", "d", "Ø" + tall(2 * r)), ("Toppsirkel", "dₐ", "Ø" + tall(2 * ra)),
           ("Rotsirkel", "df", "Ø" + tall(2 * rf)), ("Tannbredde", "b", tall(b)),
           ("Tanntykkelse, slakk", "", f"–{tall(p['slakk'])}"),
           ("Samvirkende hjul", "", "z = 73, KERS-TA-1"), ("Akseavstand", "a", "36")]
    rad = [(x, y.replace("d_f", "d_f"), z) for x, y, z in rad]
    ark.tabell(212, 120, rad, (35, 8, 32), h=5.0, tittel="Tannhjulsdata")
    merk = ["Tennene er vist forenklet (ISO 2203).", "Evolventtenner, rotfoten radiell ned til rotsirkelen.",
            "Printes med navenden ned, tennene opp.",
            "Klemnavet strammes med 2 × M3 × 12 (ISO 7045)", "og mutter M3 (ISO 4032)."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Pinjong z17, m 0,8", "KERS-PI-1", "PLA", masse(ekte), "2:1")
    ark.lagre(UT / "KERS-PI-1.pdf")



# ==========================================================================================
# Ramme RA-1 (generativ design, STL fra Fusion)
# ==========================================================================================

def tegn_ramme():
    import KERS_GD_oppsett as gd
    m = tv.les_stl(ROT / "Bygg" / "Ramme v1.stl")
    rb, rl = gd.D_BORING / 2, gd.D_LEPPE / 2
    xa = (gd.X_BAK, gd.X_SV, gd.X_FOR)
    za = gd.Z_A
    h0, h1 = gd.Y_HUS
    k = 1.0
    ark = Ark()
    kw = dict(vinkel=35, tol=0.03)
    # Sideriss sett fra -y (trommelsiden)
    S = Riss(ark, m, X, Z, (108, 126), k, **kw)
    sp = lambda x, z: S.p(x, -h1, z)
    for x in xa:
        ark.kryss(sp(x, za), k * rb, ut=3)
    yb = sp(0, 47.2)[1] + 4
    ark.maal_h(sp(xa[0], za + rb + 3), sp(xa[1], za + 6), yb, tall(xa[1] - xa[0]))
    ark.maal_h(sp(xa[0], za + rb + 3), sp(xa[2], za + rb + 3), yb + 6, tall(xa[2] - xa[0]))
    ark.maal_v(sp(gd.BJELKE_BAK[0], gd.BJELKE_BAK[2]), sp(xa[0] - rb - 3, za), sp(0, 0)[0] - 2, tall(za - gd.BJELKE_BAK[2]))
    ark.henvisning(sp(xa[2] + rb * 0.7, za + rb * 0.7), sp(xa[2] + 14, za + 18), f"4 × Ø{tall(2 * rb)}")
    ark.henvisning(sp(xa[1] - rb * 0.7, za - rb * 0.7), sp(xa[1] - 18, za - 22), f"U-sete Ø{tall(2 * rb)}", ha="right")
    ark.maal_v(sp(gd.U_X[1], za), sp(gd.U_X[1], gd.U_Z[1]), sp(gd.U_X[1], 0)[0] + 4, tall(gd.U_Z[1] - za))
    ark.snittplan(sp(xa[0], 8), sp(xa[0], 50), [-1, 0], "A")
    # Riss ovenfra
    O = Riss(ark, m, X, Y, (108, 84), k, **kw)
    op = lambda x, y: O.p(x, y, 47.2)
    for x in xa:
        ark.senterlinje(op(x, -h1 - 2), op(x, h1 + 2))
    for xs in gd.X_SKRUER:
        for ys in (gd.Y_SKRUE, -gd.Y_SKRUE):
            ark.kryss(op(xs, ys), k * 1.4, ut=1.5)
    ark.maal_v(op(xa[0], -h1), op(xa[0], h1), op(0, 0)[0] - 4, tall(2 * h1))
    ark.maal_h(op(gd.X_SKRUER[0], gd.Y_SKRUE), op(gd.X_SKRUER[1], gd.Y_SKRUE), op(0, h1)[1] + 4,
               tall(gd.X_SKRUER[1] - gd.X_SKRUER[0]))
    ark.henvisning(op(gd.X_SKRUER[1] + 0.9, gd.Y_SKRUE + 0.9), op(gd.X_SKRUER[1] + 24, h1 + 6),
                   f"4 × Ø{tall(gd.M3_KJERNE)}, dybde 10")
    # Snitt A-A gjennom bakakselen (x = 34), sett mot -x, til venstre for sideriss
    A = Riss(ark, m, Y, Z, (60, 126), k, snitt=(xa[0], 0, 0), **kw)
    a = lambda y, z: A.p(xa[0], y, z)
    ark.overskrift(a(0, 0)[0], sp(0, 47.2)[1] + 12, "A–A")
    ark.senterlinje(a(-h1 - 3, za), a(h1 + 3, za))
    ark.maal_h(a(-h1, za + rb), a(h1, za + rb), a(0, za + 14)[1], tall(2 * h1))
    ark.maal_h(a(-h0, za), a(h0, za), a(0, za - 18)[1], tall(2 * h0), )
    ark.maal_v(a(h1, za - rb), a(h1, za + rb), a(h1, 0)[0] + 4, "Ø" + tall(2 * rb))
    ark.henvisning(a(-h0 - 0.5, za - rl), a(-h0 + 10, za - 22), f"Leppe Ø{tall(2 * rl)} × 1")
    merk = ["Fri form fra generativ design (Fusion).",
            "Formen er gitt av STL-filen «Ramme v1.stl».",
            "Målsatte flater er bevart geometri.",
            "Printes liggende på siden, 20 % fyll, support.",
            "Lagerboringer Ø19,1 etter prøveringene."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Ramme", "KERS-RA-1", "PLA", "23 g (20 % fyll)", "1:1")
    ark.lagre(UT / "KERS-RA-1.pdf")



FUSION = ROT / "Bygg" / "Fusion_eksport"  # STL eksportert fra Fusion (Test-stasjon), 2026-11-17


def hjulaksel():
    """Drivakselen fra Fusion, flyttet så aksen ligger langs z med midten i origo."""
    return tv.les_stl(FUSION / "Drivaksling.stl").rotate([-90, 0, 0]).translate([-34, -37, 42.5])


def tegn_aksel():
    """Hjulaksel (printet): Ø6 mellom lagrene og firkant i hver ende for hjulene."""
    m = hjulaksel()
    zf, lf = 32.5, 10.0
    s2 = float(np.asarray(m.slice(42.0).to_polygons()[0])[:, 0].max())  # halv firkantbredde
    k = 2.0
    ark = Ark()
    F = Riss(ark, m, Z, X, (152, 135), k)
    f = lambda z, x=0.0: F.p(x, 0, z)
    F.senterlinje((0, 0, -zf - lf - 3), (0, 0, zf + lf + 3))
    ark.maal_h(f(-zf - lf, s2), f(zf + lf, s2), f(0, 3)[1] + 14, tall(2 * (zf + lf)))
    ark.maal_h(f(-zf, 3), f(zf, 3), f(0, 3)[1] + 7, tall(2 * zf))
    ark.maal_h(f(zf, -s2), f(zf + lf, -s2), f(0, -3)[1] - 7, tall(lf))
    ark.maal_v(f(0, -3), f(0, 3), f(-12, 0)[0], "Ø" + tall(6))
    E = Riss(ark, m, -Y, X, (40, 135), k)
    e0 = E.p(0, 0, 0)
    ark.kryss(e0, k * 3)
    ark.maal_h(E.p(-s2, s2, 0), E.p(s2, -s2, 0), e0[1] + k * 3 + 6, tall(2 * s2))
    ark.maal_v(E.p(-s2, -s2, 0), E.p(s2, -s2, 0), e0[0] + k * 3 + 6, tall(2 * s2))
    merk = ["2 stk. (bakaksel og foraksel).",
            "Firkanten passer firkanthullet 4,3 i felgen",
            "(KERS-FE-1). Ø6 går gjennom innerringen i",
            "626-lagrene. Modell fra Fusion (Drivaksling)."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Hjulaksel", "KERS-AK-1", "PLA", masse(m), "2:1")
    ark.lagre(UT / "KERS-AK-1.pdf")


def tannhjul73():
    """Tannhjul z73 forenklet (toppsylinder i stedet for tenner), mål fra Fusion-modellen:
    tenner z 0-5, nav Ø14 z 5-11, boring Ø6,2 og radielt hull Ø2,5 mot -y midt i navet."""
    hull = mf.Manifold.cylinder(6, 1.25, 1.25, 48).rotate([90, 0, 0]).translate([0, -2, 8])
    return T.syl(30, 0, 5) + T.syl(7, 5, 11) - T.syl(3.1, -1, 12) - hull


def tegn_tannhjul():
    ekte = tv.les_stl(FUSION / "Spur_Gear_73_teeth.stl")
    m = tannhjul73()
    z, mod = 73, 0.8
    r, ra = mod * z / 2, 30.0
    rf = ekte.slice(2.5).to_polygons()
    rf = min(np.linalg.norm(np.asarray(p), axis=1).min() for p in rf if len(p) > 100)
    k = 1.0
    ark = Ark()
    # Riss fra tannsiden (z = 0), hullet opp
    F = Riss(ark, m, X, -Y, (72, 120), k)
    f0 = F.p(0, 0, 0)
    ark.kryss(f0, k * ra)
    ark.sirkel(f0, k * r, lw=tv.TYNN, ls=tv.STREKPUNKT)
    ark.sirkel(f0, k * rf, lw=tv.TYNN)
    ark.snittlinje(f0 + [0, k * ra + 3], [1, 0], "A", [0, 1])
    ark.snittlinje(f0 - [0, k * ra + 3], [1, 0], "A", [0, -1])
    # Snitt A-A gjennom hullet, tannsiden til høyre
    A = Riss(ark, m, -Z, -Y, (128, 120), k, snitt=(0, 0, 0))
    a = lambda zz, rr: A.p(0, -rr, zz)
    A.senterlinje((0, 0, 14), (0, 0, -3))
    for sg in (1, -1):
        ark.linje([a(-2, sg * r), a(7, sg * r)], lw=tv.TYNN, ls=tv.STREKPUNKT, z=5)
    ark.overskrift(a(5.5, 0)[0], a(0, ra)[1] + 14, "A–A")
    ark.maal_h(a(0, ra), a(5, ra), a(0, ra)[1] + 4, "5")
    ark.maal_h(a(0, ra), a(11, 7), a(0, ra)[1] + 9, "11")
    ark.maal_v(a(0, -ra), a(0, ra), a(-5, 0)[0], "Ø60")
    ark.detaljmerke(a(8, 5), 8, "B", vinkel=150)
    # Detalj B (5:1): navet med hullet
    kd = 5.0
    cB_m = (0, -4.6, 7.5)
    cB = np.array([196.0, 150.0])
    D = Riss(ark, m, -Z, -Y, cB - kd * (tv.basis(-Z, -Y)[:2] @ np.asarray(cB_m)), kd, snitt=(0, 0, 0), tegn=False)
    D.tegn(klipp=(cB, kd * 5.2))
    ark.overskrift(cB[0], cB[1] + kd * 5.2 + 3, "B (5:1)")
    d = lambda zz, rr: D.p(0, -rr, zz)
    ark.senterlinje(d(8, 2.2), d(8, 8.0))
    ark.maal_h(d(8, 7.0), d(11, 7.0), d(0, 7.0)[1] + 4, "3")
    ark.henvisning(d(8.9, 6.0), cB + [16, 22], "Ø2,5")
    ark.henvisning(d(10.6, 3.1), cB + [16, -22], "Ø6,2")
    ark.henvisning(d(10.2, 7.0), cB + [-16, 24], "Ø14", ha="right")
    rad = [("Antall tenner", "z", "73"), ("Modul", "m", "0,8"), ("Trykkvinkel", "α", "20°"),
           ("Profilforskyvning", "x", "0"), ("Delesirkel", "d", "Ø" + tall(2 * r)),
           ("Toppsirkel", "dₐ", "Ø60"), ("Rotsirkel", "df", "Ø" + tall(round(2 * rf, 2))),
           ("Tannbredde", "b", "5"), ("Samvirkende hjul", "", "z = 17, KERS-PI-1"),
           ("Akseavstand", "a", "36")]
    ark.tabell(212, 112, rad, (35, 8, 32), h=5.0, tittel="Tannhjulsdata")
    merk = ["Tennene er vist forenklet (ISO 2203).",
            "Modell: Spur Gear 73 teeth i Fusion.",
            "Ø2,5-hullet er for en M3-settskrue mot akselen,",
            "skrudd rett i PLA. Printes med tannsiden ned."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Tannhjul z73, m 0,8", "KERS-TA-1", "PLA", masse(ekte), "1:1 (B 5:1)")
    ark.lagre(UT / "KERS-TA-1.pdf")


def magnetholder():
    """Magnetholderen fra Fusion: aksen langs y flyttes til z, lommesiden på z = 0."""
    return tv.les_stl(FUSION / "Magnetholder.stl").rotate([-90, 0, 0]).translate([0, 0, 29.5]).mirror([0, 0, 1]).translate([0, 0, 3])


def tegn_magnetholder():
    m = magnetholder()
    k = 2.0
    ark = Ark()
    rl, cl, dl = 1.6, 9.05, 2.2
    # Riss mot lommesiden (z = 0), betrakteren på -z
    F = Riss(ark, m, X, -Y, (80, 120), k)
    f0 = F.p(0, 0, 0)
    ark.kryss(f0, k * 15)
    for sg in (1, -1):
        ark.kryss(F.p(0, sg * cl, 0), k * rl, ut=2)
    pa, pb = F.p(0, -cl, 0), F.p(0, cl, 0)
    ark.maal_v(pa, pb, f0[0] - k * 15 - 8, tall(2 * cl))
    ark.maal_dia(f0, k * 15, -35, "Ø30")
    ark.henvisning(F.p(rl * 0.7, -cl - rl * 0.7, 0), f0 + [26, 42], f"2 × Ø{tall(2 * rl)}, dybde {tall(dl)}")
    ark.snittlinje(f0 + [0, k * 15 + 3], [1, 0], "A", [0, 1])
    ark.snittlinje(f0 - [0, k * 15 + 3], [1, 0], "A", [0, -1])
    # Snitt A-A gjennom lommene, lommesiden til høyre
    A = Riss(ark, m, -Z, -Y, (160, 120), k, snitt=(0, 0, 0))
    a = lambda zz, rr: A.p(0, -rr, zz)
    A.senterlinje((0, 0, 6), (0, 0, -3))
    ark.overskrift(a(1.5, 0)[0], a(0, 15)[1] + 13, "A–A")
    ark.maal_h(a(0, 15), a(3, 15), a(0, 15)[1] + 4, "3")
    ark.maal_v(a(3, -3.1), a(3, 3.1), a(3, 0)[0] - 6, "Ø6,2")
    ark.henvisning(a(1.6, 3.1 + 3 - (9 - 1.4 ** 2) ** 0.5), a(5, 9) + [-6, 8], "R3", ha="right")
    merk = ["Klemmes mellom svinghjulet og pinjongnavet.",
            "R3-siden vender mot svinghjulet (sveisen).",
            "To magneter Ø3 × 2 limes i lommene",
            "(Hall-sensor). Modell: Magnetholder i Fusion."]
    ark.merknader(24, 12 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.tittelfelt("Magnetholder", "KERS-MH-1", "PLA", masse(m), "2:1")
    ark.lagre(UT / "KERS-MH-1.pdf")


# ==========================================================================================
# Sammenstillinger
# ==========================================================================================

def plasser(m, x0, y0, z0, s, dreid=0.0):
    """Del med aksen langs lokal z settes inn med aksen langs bilens y. s = +1: lokal z mot +y,
    s = -1: mot -y. Lokal x blir bilens x. dreid: rotasjon om egen akse først (grader)."""
    if dreid:
        m = m.rotate([0, 0, dreid])
    return m.transform(np.array([[1, 0, 0, x0], [0, 0, s, y0], [0, -s, 0, z0]], float))


def bil_deler(trommel="TR-0", hylse_pluss=6.0, dreid_nav=90):
    """Alle delene i bilens koordinater (x fra bakkant, y fra midtlinjen, +y tannhjulssiden, z fra
    bakken). Tannhjulene vises forenklet (toppsylinder), slik ISO 2203 gjør i riss."""
    import KERS_GD_oppsett as gd
    xa, za = (gd.X_BAK, gd.X_SV, gd.X_FOR), gd.Z_A
    d = {}
    d["ramme"] = tv.les_stl(ROT / "Bygg" / "Ramme v1.stl")
    lokk = tv.les_stl(ROT / "Bygg" / "Lagerlokk -y.stl")
    d["lokk"] = [lokk, lokk.mirror([0, 1, 0])]
    lager = T.syl(9.5, 0, 6) - T.syl(3.0, -1, 7)
    d["lager"] = [plasser(lager, x, s * 26.5, za, s) for x in xa for s in (1, -1)]
    d["svinghjul"] = plasser(T.syl(30, -3, 3) - T.syl(3, -4, 4), 70, 0, za, 1)
    d["aksel_sv"] = plasser(T.syl(3, -32.5, 32.5), 70, 0, za, 1)
    p0 = T.P_TR0 if trommel == "TR-0" else T.P
    d["trommel"] = plasser(T.trommel(p0, 8.0 if trommel == "TR-0" else 12.0), 70, -3, za, -1, dreid=dreid_nav)
    pin, _ = pinjong_forenklet(PI.P)
    # magnetholderen ligger mot svinghjulet (R3-siden mot sveisen), pinjongen utenfor den
    d["pinjong"] = plasser(pin, 70, 6, za, 1, dreid=dreid_nav)
    dreid_mag = 90 - dreid_nav  # lommene i snittplanet
    d["magnetholder"] = plasser(magnetholder(), 70, 6, za, -1, dreid=dreid_mag)
    magneter = T.syl(1.5, 0, 2).translate([0, 9.05, 0]) + T.syl(1.5, 0, 2).translate([0, -9.05, 0])
    d["magneter"] = plasser(magneter, 70, 6, za, -1, dreid=dreid_mag)
    d["tannhjul"] = plasser(tannhjul73(), gd.X_BAK, 5, za, 1)
    d["aksler"] = [plasser(hjulaksel(), x, 0, za, 1) for x in (gd.X_BAK, gd.X_FOR)]
    p = H.P
    felg = H.hjul(p, True)
    r_sete = p["d_ytre"] / 2 - p["t_dekk"]
    ri = r_sete - p["strekk"] / 2
    zd = (p["bredde"] - p["b_dekk"]) / 2 + 0.1
    dekk = H.dekk(p).scale([r_sete / ri, r_sete / ri, 1]).translate([0, 0, zd])
    d["felg"] = [plasser(felg, x, s * 41.5, za, -s) for x in (gd.X_BAK, gd.X_FOR) for s in (1, -1)]
    d["dekk"] = [plasser(dekk, x, s * 41.5, za, -s) for x in (gd.X_BAK, gd.X_FOR) for s in (1, -1)]
    lt = 20.4 if trommel == "TR-0" else 20.5
    d["hylse"] = [plasser(AH.hylse(AH.P, 3.0), 70, -3 - lt - 0.1, za, -1),
                  plasser(AH.hylse(AH.P, hylse_pluss), 70, 26.5 - hylse_pluss, za, 1)]
    skrue = T.syl(2.8, 0, 1.8) + T.syl(1.5, -16, 0)
    d["skruer_lokk"] = [skrue.translate([xs, s * gd.Y_SKRUE, gd.U_Z[1] + gd.H_LOKK]) for xs in gd.X_SKRUER
                        for s in (1, -1)]
    return d


STYKKLISTE = [
    (1, 1, "Ramme", "PLA", "KERS-RA-1"),
    (2, 2, "Lagerlokk", "PLA", "KERS-LL-1"),
    (3, 1, "Svinghjul med aksel", "Stål, sveist", "Utdelt"),
    (4, 1, "Snortrommel TR-0", "PLA", "KERS-TR-0"),
    (5, 1, "Pinjong z17", "PLA", "KERS-PI-1"),
    (6, 1, "Magnetholder", "PLA", "KERS-MH-1"),
    (7, 1, "Tannhjul z73", "PLA", "KERS-TA-1"),
    (8, 2, "Hjulaksel", "PLA", "KERS-AK-1"),
    (9, 4, "Felg", "PLA", "KERS-FE-1"),
    (10, 4, "Dekk", "TPU", "KERS-DE-1"),
    (11, 2, "Avstandshylse", "PLA", "KERS-AH-1, L = 3,0 og 6,0"),
    (12, 6, "Kulelager 626", "6 × 19 × 6", "Innkjøpt"),
    (13, 4, "Linsehodeskrue M3 × 12", "ISO 7045", "Klemnav, pos. 4 og 5"),
    (14, 4, "Mutter M3", "ISO 4032", "Klemnav, pos. 4 og 5"),
    (15, 4, "Linsehodeskrue M3 × 16", "ISO 7045", "Lagerlokk, pos. 2"),
    (16, 1, "Settskrue M3", "ISO 4027", "I navet på pos. 7"),
    (17, 2, "Magnet, neodym", "Ø3 × 2", "Limt i pos. 6"),
]


def tegn_sm1():
    import KERS_GD_oppsett as gd
    d = bil_deler()
    za = gd.Z_A
    deler = [dict(m=d["ramme"], skravur="45", tetthet=3, vinkel=35)]
    deler += [dict(m=m, skravur="-45", tetthet=4) for m in d["lokk"]]
    deler += [dict(m=m, skravur="lager") for m in d["lager"]]
    deler += [dict(m=d["svinghjul"], skravur="-45", tetthet=6), dict(m=d["aksel_sv"], skravur=None),
              dict(m=d["trommel"], skravur="45", tetthet=6), dict(m=d["pinjong"], skravur="45", tetthet=6),
              dict(m=d["tannhjul"], skravur="45"), dict(m=d["magnetholder"], skravur="-45", tetthet=4),
              dict(m=d["magneter"], skravur="45", tetthet=8)]
    deler += [dict(m=m, skravur=None) for m in d["aksler"]]
    deler += [dict(m=m, skravur="-45") for m in d["felg"]] + [dict(m=m, skravur="45", tetthet=8) for m in d["dekk"]]
    deler += [dict(m=m, skravur="-45", tetthet=6) for m in d["hylse"]]
    deler += [dict(m=m, skravur=None) for m in d["skruer_lokk"]]
    k = 1.0
    ark = Ark(liggende=False)
    # Sideriss sett fra tannhjulssiden (+y): fronten til venstre
    S = tv.Sammenstilling(ark, deler, -X, Z, (110, 243), k, ref=(70, 0, za), tol=0.03)
    for x in (gd.X_BAK, gd.X_SV, gd.X_FOR):
        ark.kryss(S.p(x, 50, za), 9, ut=2)
    ark.snittplan(S.p(70, 50, 1), S.p(70, 50, 53), [1, 0], "A")
    # Snitt A-A gjennom svinghjulsakselen, sett bakover (-x)
    A = tv.Sammenstilling(ark, deler, Y, Z, (110, 155), k, ref=(70, 0, za), snitt=(70, 0, 0), tol=0.03)
    ark.overskrift(92, 192, "A–A")
    ark.senterlinje(A.p(70, -46, za), A.p(70, 46, za))
    # posisjonsnumre: (riss, modellpunkt, ballongens sentrum på arket, pos)
    xs1 = gd.X_SKRUER[0]
    for R, q, c, nr in (
            (S, (124, 0, 44), (34, 280), 1), (S, (70, 32.5, 46.5), (96, 282), 2), (S, (70, 3, 6.5), (140, 207), 3),
            (S, (34, 43, 35.2), (190, 263), 8), (S, (121, 41.5, 33), (30, 229), 9), (S, (34, 41, 67.6), (176, 282), 10),
            (A, (70, -16, za + 9.4), (40, 189), 4), (A, (70, 10, za + 7.2), (152, 191), 5),
            (A, (70, 4.5, za - 12.5), (138, 133), 6), (A, (34, 7.5, 9), (162, 125), 7),
            (A, (70, -25, za + 5.5), (40, 173), 11),
            (A, (70, 29.5, za + 9.0), (188, 167), 12), (A, (70, -20.5, za - 7), (40, 147), 13),
            (A, (70, -18.5, za - 9.6), (40, 131), 14), (A, (xs1, 29.5, 48.6), (184, 189), 15),
            (A, (70, -29.0, 46.0), (60, 196), 2), (A, (70, 1.5, za + 18), (128, 196), 3)):
        ark.ballong(R.p(q), c, nr)
    ark.ballong(None, (170, 125), 16)  # settskruen i tannhjulsnavet (skjult), gruppert med pos. 7
    ark.ballong(None, (146, 133), 17)  # magnetene i holderen, gruppert med pos. 6
    ark.tittelfelt("KERS-bil, sammenstilling", "KERS-SM-1", "Se stykkliste", "270 g", "1:1",
                   type_="Sammenstillingstegning")
    ark.stykkliste([tuple(map(str, r)) for r in STYKKLISTE], bredder=(10, 10, 52, 34, 74), h=4.2)
    ark.merknader(158, 154, ["Tannhjulene er vist forenklet", "(ISO 2203), lagrene forenklet", "(ISO 8826-1)."], h=2.2)
    ark.lagre(UT / "KERS-SM-1.pdf")




def tegn_sm2():
    """Svinghjulsakselen som ferdig montert enhet (tres på før den senkes ned i U-setene)."""
    import KERS_GD_oppsett as gd
    d = bil_deler(dreid_nav=0)
    za = gd.Z_A
    lag = [m for m in d["lager"] if abs(m.bounding_box()[0] - 60.5) < 1]
    deler = [dict(m=d["svinghjul"], skravur="-45", tetthet=6), dict(m=d["aksel_sv"], skravur=None),
             dict(m=d["trommel"], skravur="45", tetthet=4), dict(m=d["pinjong"], skravur="45", tetthet=4),
             dict(m=d["magnetholder"], skravur="-45", tetthet=4), dict(m=d["magneter"], skravur="45", tetthet=8)]
    deler += [dict(m=m, skravur="-45", tetthet=4) for m in d["hylse"]] + [dict(m=m, skravur="lager") for m in lag]
    k = 2.0
    ark = Ark(liggende=False)
    # Snitt A-A i planet gjennom akselen (z = 34), sett ovenfra; +y (tannhjulssiden) til høyre
    A = tv.Sammenstilling(ark, deler, Y, -X, (110, 186), k, ref=(70, 0, za), snitt=(0, 0, za))
    p = lambda y, x=0.0: A.p(70 + x, y, za)
    ark.senterlinje(p(-36), p(36))
    # monteringsmål langs akselen, fra svinghjulets midtplan
    yt = p(0, -30)[1] + 6
    ark.maal_h(p(-32.5, -3), p(32.5, -3), yt + 12, tall(65))
    ark.maal_h(p(-26.5, -9.5), p(26.5, -9.5), yt + 6, tall(53))
    yb = p(0, 30)[1] - 6
    ark.maal_h(p(0, 3), p(6, 7.6), yb, tall(6))
    ark.maal_h(p(0, 3), p(-3, 9), yb, tall(3))
    for q, c, nr in (((0.0, -22), (70, 262), 1), ((-12, -8.5), (40, 248), 2), ((12, -7.4), (150, 254), 3),
                     ((4.5, -13), (96, 280), 4), ((5.0, -9.05), (124, 280), 5),
                     ((-25.0, -5.5), (32, 220), 6), ((24.0, -5.5), (188, 222), 7), ((29.5, -7.5), (188, 240), 8),
                     ((-29.5, -7.5), (32, 236), 8), ((-20, 6.0), (40, 150), 9), ((-20, 7.5), (60, 140), 10)):
        ark.ballong(A.p(70 + q[1], q[0], za), c, nr)
    rader = [(1, 1, "Svinghjul med aksel", "Stål, sveist", "Utdelt (pos. 3 i SM-1)"),
             (2, 1, "Snortrommel TR-0", "PLA", "KERS-TR-0"),
             (3, 1, "Pinjong z17", "PLA", "KERS-PI-1"),
             (4, 1, "Magnetholder", "PLA", "KERS-MH-1"),
             (5, 2, "Magnet, neodym", "Ø3 × 2", "Limt i pos. 4"),
             (6, 1, "Avstandshylse L = 3,0", "PLA", "KERS-AH-1"),
             (7, 1, "Avstandshylse L = 6,0", "PLA", "KERS-AH-1"),
             (8, 2, "Kulelager 626", "6 × 19 × 6", "Innkjøpt"),
             (9, 4, "Linsehodeskrue M3 × 12", "ISO 7045", "I klemnavene, pos. 2 og 3"),
             (10, 4, "Mutter M3", "ISO 4032", "I klemnavene, pos. 2 og 3")]
    ark.tittelfelt("Svinghjulsaksel, montert", "KERS-SM-2", "Se stykkliste", "–", "2:1",
                   type_="Sammenstillingstegning")
    yt_ = ark.stykkliste([tuple(map(str, r)) for r in rader], bredder=(10, 10, 52, 34, 74), h=5.0)
    merk = ["Snitt i planet gjennom akselen, sett ovenfra. Mål fra svinghjulets midtplan.",
            "Delene tres på fra akselendene og klemmes med navene, så senkes",
            "enheten ned i U-setene i rammen og lagerlokkene skrus på (SM-1).",
            "Hylsene velges så akselen får 0,1–0,3 mm klaring sideveis."]
    ark.merknader(24, yt_ + 4 + 4.2 * (len(merk) - 1), merk, h=2.3)
    ark.lagre(UT / "KERS-SM-2.pdf")


TEGNINGER = {
    "TR-0": lambda: tegn_trommel("TR-0"),
    "TR-1": lambda: tegn_trommel("TR-1"),
    "FE-1": tegn_felg,
    "DE-1": tegn_dekk,
    "AH-1": tegn_hylse,
    "LL-1": tegn_lokk,
    "PI-1": tegn_pinjong,
    "RA-1": tegn_ramme,
    "AK-1": tegn_aksel,
    "TA-1": tegn_tannhjul,
    "MH-1": tegn_magnetholder,
    "SM-1": tegn_sm1,
    "SM-2": tegn_sm2,
}


def main():
    UT.mkdir(parents=True, exist_ok=True)
    valg = sys.argv[1:] or list(TEGNINGER)
    for nr in valg:
        TEGNINGER[nr]()
        print("Skrevet", UT / f"KERS-{nr}.pdf")


if __name__ == "__main__":
    main()
