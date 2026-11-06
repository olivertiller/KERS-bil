"""Friksjonshjul til drillen, for å snurre opp svinghjulet før utsnurringstesten (3MF i mm).

    python3 Code/drillhjul.py                    # standardstørrelsene
    python3 Code/drillhjul.py --diametre 80 120  # andre diametre

Hjulet spennes i chucken med en sekskantstamme og trykkes ned mot toppen av svinghjulet. Kanten
av svinghjulet stikker 14 mm over veggene (topp z = 64, vegg z = 50). Drivhjulet går mellom y = ±5,
klar av pinjong (y 5 til 19,5) og trommel (y -3 og ned, men den sitter lavt). Utveksling:

    n_svinghjul = n_drill * D_drivhjul / 60

Makita DDF483 uten last: trinn 1 = 500 o/min, trinn 2 = 1700 o/min (litt lavere under last).

Skriver til Bygg/Drillhjul/, for hver diameter D:
  drillhjul_D<D>_PLA.3mf   nav, eiker, felg og sekskantstamme i ett. Printes med eikene ned, stammen opp.
  drillhjul_D<D>_dekk_TPU.3mf   dekk som strekkes på felgen
  drillhjul_forhandsvisning.png
"""

import argparse
import math
from pathlib import Path

import manifold3d as mf
import numpy as np

from hjul import eikeskive
from snortrommel import RHO_PLA, skriv_3mf, syl

UT = Path(__file__).resolve().parent.parent / "Bygg" / "Drillhjul"
RHO_TPU = 1.21e-3
D_SVINGHJUL = 60.0
DIAMETRE = [60.0, 100.0, 140.0]
N_DRILL = {"trinn 1": 500, "trinn 2": 1700}

P = dict(
    bredde=10.0,  # felg (svinghjulet er 6 bredt)
    t_dekk=2.0,
    b_dekk=8.0,
    strekk_rel=0.01,  # dekkets indre Ø er 1 % mindre enn setet
    skulder=1.0,
    t_felg=3.0,
    t_skive=3.0,
    n_eiker=6,
    b_eike=6.0,
    d_nav=20.0,
    sekskant=10.0,  # nøkkelvidde for stammen (chucken tar 1,5 til 13)
    l_stamme=25.0,
    d_rot=16.0,  # avrundet overgang mellom nav og stamme
    l_rot=4.0,
    fas=0.8,
)


def sekskant(nv, z0, z1, fas=0.0):
    """Sekskantprisme med nøkkelvidde nv langs z, med fas i toppen."""
    r = nv / math.sqrt(3)  # omskrevet radius
    m = mf.Manifold.cylinder(z1 - z0 - fas, r, r, 6).translate([0, 0, z0])
    if fas:
        m = m + mf.Manifold.cylinder(fas, r, r - fas, 6).translate([0, 0, z1 - fas])
    return m


def drillhjul(p, d):
    b, kant = p["bredde"], (p["bredde"] - p["b_dekk"]) / 2
    r_sete = d / 2 - p["t_dekk"]
    rs = r_sete + p["skulder"]
    r_inn = r_sete - p["t_felg"]
    felg = syl(r_sete, 0, b) + syl(rs, 0, kant) + syl(rs, b - kant, b) - syl(r_inn, -1, b + 1)
    rn = p["d_nav"] / 2
    m = felg + eikeskive(p, rn - 0.5, r_inn + 0.5) + syl(rn, 0, b)
    # overgang nav -> stamme (kjegle) og sekskantstammen
    rr = p["d_rot"] / 2
    m = m + syl(rr, b, b + p["l_rot"], p["sekskant"] / 2)
    m = m + sekskant(p["sekskant"], b, b + p["l_rot"] + p["l_stamme"], p["fas"])
    return m


def dekk(p, d):
    r_sete = d / 2 - p["t_dekk"]
    ri = r_sete * (1 - p["strekk_rel"])
    ro = ri + p["t_dekk"]
    bd, f = p["b_dekk"] - 0.2, 0.5
    m = syl(ro - f, 0, f, ro) + syl(ro, f, bd - f) + syl(ro, bd - f, bd, ro - f)
    return m - syl(ri, -1, bd + 1)


def tegn(p, deler, sti):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    n = len(deler)
    fig, ax = plt.subplots(1, n, figsize=(4.6 * n, 5.4))
    ax = np.atleast_1d(ax)
    for a, (d, (hj, dk)) in zip(ax, deler.items()):
        kant = (p["bredde"] - p["b_dekk"] + 0.2) / 2
        for m, fc in ((hj, "#d9d9d9"), (dk.translate([0, 0, kant]), "#2b2b2b")):
            for poly in m.rotate([90, 0, 0]).slice(0.0).to_polygons():
                pts = np.asarray(poly) * [1, -1]
                a.add_patch(Polygon(pts[:, [1, 0]], fc=fc, ec="k", lw=0.5))
        # svinghjulet under, i kontakt
        r = d / 2
        a.add_patch(Polygon([[2, -r], [8, -r], [8, -r - 8], [2, -r - 8]], fc="#8a8f99", ec="k", lw=0.5))
        a.text(9, -r - 5, "svinghjul Ø60", fontsize=7)
        a.set_xlim(-6, 45); a.set_ylim(-r - 10, r + 3); a.set_aspect("equal")
        tekst = ", ".join(f"{k}: {v * d / D_SVINGHJUL:.0f}" for k, v in N_DRILL.items())
        a.set_title(f"Ø{d:g}  (svinghjul o/min, {tekst})", fontsize=9)
        a.set_xlabel("z [mm] (0 = eikesiden, ned på printeren)")
    fig.tight_layout()
    fig.savefig(sti, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diametre", type=float, nargs="+", default=DIAMETRE, help="ferdig ytre Ø med dekk [mm]")
    a = ap.parse_args()
    UT.mkdir(parents=True, exist_ok=True)

    deler = {}
    print(f"{'Ø':>5}  {'PLA':>6}  {'TPU':>5}   svinghjulets turtall (drill uten last)")
    for d in a.diametre:
        hj, dk = drillhjul(P, d), dekk(P, d)
        for m, navn in ((hj, f"drillhjul_D{d:g}_PLA"), (dk, f"drillhjul_D{d:g}_dekk_TPU")):
            assert m.status() == mf.Error.NoError, navn
            skriv_3mf(m, UT / f"{navn}.3mf", navn)
        deler[d] = (hj, dk)
        turtall = ", ".join(f"{k} {v * d / D_SVINGHJUL:5.0f}" for k, v in N_DRILL.items())
        print(f"{d:5.0f}  {hj.volume() * RHO_PLA:5.1f}g  {dk.volume() * RHO_TPU:4.1f}g   {turtall} o/min")
    tegn(P, deler, UT / "drillhjul_forhandsvisning.png")
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
