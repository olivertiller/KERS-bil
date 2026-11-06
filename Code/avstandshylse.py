"""Avstandshylser som holder svinghjulsakselen på plass sideveis (3MF i mm).

    python3 Code/avstandshylse.py                    # hele settet, 1,0 til 7,0 mm
    python3 Code/avstandshylse.py --lengder 2.8 4.6  # bare de du har målt

Én hylse på hver side av svinghjulet:
  -y: mellom snortrommelens klemnav og lageret (nominelt 3,0 mm luft)
  +y: mellom magnetholderen/pinjongnavet og lageret (nominelt 7,0 mm fra pinjongnavet, minus
      magnetholderen)
Begge sider bruker samme hylse, bare lengden er ulik. Mål luften med skyvelæret (dybdepinnen) med
svinghjulet sentrert, og velg den nærmeste hylsen som er 0,1 til 0,2 mm kortere. Trommel og pinjong
har klemnav, så resten kan justeres ved å skyve navet inntil hylsen før du strammer.

Hylsen er C-formet, så den klikkes på akselen uten å ta den ut av lagrene. Åpningen (5,3 mm) er
smalere enn akselen (Ø6), og skråkanter leder den inn. Kragen Ø7,8 ligger an mot lagerets innerring
(som roterer med akselen) og aldri mot tetningen eller ytterringen. Den brede delen Ø12 ligger an mot
navet. Printes stående med den brede enden ned og kragen opp, uten support.

Skriver til Bygg/Avstandshylser/: hylse_L<lengde>.3mf og hylse_forhandsvisning.png
"""

import argparse
import math
from pathlib import Path

import manifold3d as mf
import numpy as np

from snortrommel import skriv_3mf

UT = Path(__file__).resolve().parent.parent / "Bygg" / "Avstandshylser"
SEG = 96

P = dict(
    boring=6.1,  # litt grep rundt Ø6-akselen
    d_ytre=12.0,  # mot navet (trommel- og pinjongnav er Ø20)
    d_krage=7.8,  # mot lagerets innerring, maks Ø8 (samme som kragen på hjulene)
    l_krage=0.8,
    apning=5.3,  # bredde på åpningen inn mot boringen (< Ø6, så den klikker på)
    apning_ytre=8.0,  # bredde ute ved Ø12, skråkanter leder akselen inn
)


def c_profil(r_ut, p):
    """C-formet ring: sirkel minus boring minus kileformet åpning mot +x."""
    ring = mf.CrossSection.circle(r_ut, SEG) - mf.CrossSection.circle(p["boring"] / 2, SEG)
    a, b = p["apning"] / 2, p["apning_ytre"] / 2
    x0 = math.sqrt((p["boring"] / 2) ** 2 - a ** 2) - 0.2  # start litt inne i boringen
    kile = mf.CrossSection([[(x0, -a), (r_ut + 1, -b), (r_ut + 1, b), (x0, a)]])
    return ring - kile


def hylse(p, lengde):
    lk = min(p["l_krage"], lengde / 3)
    kropp = c_profil(p["d_ytre"] / 2, p).extrude(lengde - lk)
    krage = c_profil(p["d_krage"] / 2, p).extrude(lk).translate([0, 0, lengde - lk])
    return kropp + krage  # printes med den brede enden ned, kragen opp


def tegn(p, deler, sti):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw=dict(height_ratios=[2, 1]))
    m = next(iter(deler.values()))
    for poly in m.slice(0.5).to_polygons():
        a0.add_patch(Polygon(np.asarray(poly), fc="#d9d9d9", ec="k", lw=0.6))
    a0.add_patch(plt.Circle((0, 0), 3.0, fc="none", ec="#b4541a", ls="--", lw=0.8))
    a0.text(0, -0.4, "aksel Ø6", ha="center", fontsize=7, color="#b4541a")
    a0.set_xlim(-7.5, 7.5); a0.set_ylim(-7.5, 7.5); a0.set_aspect("equal")
    a0.set_title("Sett ovenfra (C-form, klikkes på)")
    x = 0.0
    for lengde, m in deler.items():
        for poly in m.rotate([0, 0, 90]).rotate([-90, 0, 0]).slice(0.0).to_polygons():  # snitt på tvers av åpningen
            pts = np.asarray(poly)
            a1.add_patch(Polygon(pts + [x, 0], fc="#d9d9d9", ec="k", lw=0.4))
        a1.text(x, -1.2, f"{lengde:g}", ha="center", fontsize=8)
        x += 20
    a1.set_xlim(-8, x - 12); a1.set_ylim(-2, 8); a1.set_aspect("equal")
    a1.set_title("Snitt, lengde i mm (kragen opp, mot lageret)")
    a1.set_yticks([])
    fig.tight_layout()
    fig.savefig(sti, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lengder", type=float, nargs="+", default=list(np.arange(1.0, 7.01, 0.5)))
    a = ap.parse_args()
    UT.mkdir(parents=True, exist_ok=True)
    deler = {}
    for lengde in a.lengder:
        m = hylse(P, lengde)
        assert m.status() == mf.Error.NoError, lengde
        navn = f"hylse_L{lengde:.1f}"
        skriv_3mf(m, UT / f"{navn}.3mf", navn)
        deler[lengde] = m
    tegn(P, deler, UT / "hylse_forhandsvisning.png")
    print(f"{len(deler)} hylser ({', '.join(f'{x:g}' for x in deler)} mm) skrevet til {UT}")


if __name__ == "__main__":
    main()
