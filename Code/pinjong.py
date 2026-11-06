"""Pinjong z17, modul 0,8 med klemnav (samme feste som snortrommelen), som 3MF (mm).

    python3 Code/pinjong.py                 # standard
    python3 Code/pinjong.py --slakk 0.15    # tynnere tenner hvis inngrepet går tregt
    python3 Code/pinjong.py --boring 6.3    # hvis boringen blir for trang

Krever manifold3d (se snortrommel.py).

Pinjongen printes med navenden mot platen og tennene opp (3MF-filen er snudd slik), fordi navet
(Ø20) er større enn tannhjulet (Ø15,2) og ellers ville hengt flatt ut i luften. Koordinater i
koden: z langs aksen, z = 0 er tannhjulsflaten mot svinghjulet. I bilen er det y = 5 + z.

  * Evolventtenner, trykkvinkel 20°, uten profilforskyvning: a = 0,8 * (17 + 73) / 2 = 36 mm.
    Rotfoten er forenklet til en rett radiell linje fra grunnsirkelen ned til rotsirkelen.
    z = 17 ligger på grensen for underskjæring (17,1 ved 20°), så tennene får full profil.
  * Tannbredde 7 mm (y 5 til 12).
  * Klemnav Ø20 x 7,5 med to slisser og to M3 x 12 + mutter på tvers, 180° fra hverandre,
    akkurat som på snortrommelen (y 12 til 19,5).
"""

import argparse
import math
from pathlib import Path

import manifold3d as mf
import numpy as np

from snortrommel import P as P_TROMMEL, klemnav_kutt, skriv_3mf, syl

UT = Path(__file__).resolve().parent.parent / "Bygg" / "Pinjong"

P = dict(
    {k: P_TROMMEL[k] for k in ("boring", "d_nav", "l_nav", "slisse", "e_skrue", "m3_hull", "d_hode", "mutter", "x_sete")},
    z=17,
    modul=0.8,
    alfa=20.0,  # trykkvinkel [°]
    b=7.0,  # tannbredde (5 -> 7 2026-11-05)
    slakk=0.10,  # tannflanke-slakk: tanntykkelsen på delesirkelen reduseres med dette [mm]
)


def tannprofil(p, n_flanke=12, n_rot=6, n_topp=4):
    """Lukket polygon (mot klokka) for tannhjulssnittet."""
    z, m = p["z"], p["modul"]
    a = math.radians(p["alfa"])
    r = m * z / 2
    rg = r * math.cos(a)  # grunnsirkel
    ra = r + m  # toppsirkel
    rf = r - 1.25 * m  # rotsirkel
    inv = lambda x: math.tan(x) - x
    psi = (math.pi * m / 2 - p["slakk"]) / (2 * r)  # halv tannvinkel på delesirkelen
    theta = lambda rho: psi + inv(a) - inv(math.acos(min(1.0, rg / rho)))
    rhos = np.linspace(max(rg, rf), ra, n_flanke)
    pts = []
    for k in range(z):
        phi = 2 * math.pi * k / z
        tb = theta(max(rg, rf))
        pts.append((rf, phi - tb))
        pts += [(rho, phi - theta(rho)) for rho in rhos]
        ta = theta(ra)
        pts += [(ra, phi + t) for t in np.linspace(-ta, ta, n_topp)[1:-1]]
        pts += [(rho, phi + theta(rho)) for rho in rhos[::-1]]
        pts.append((rf, phi + tb))
        neste = 2 * math.pi * (k + 1) / z - tb
        pts += [(rf, t) for t in np.linspace(phi + tb, neste, n_rot)[1:-1]]
    return [(rho * math.cos(t), rho * math.sin(t)) for rho, t in pts], (r, rg, ra, rf)


def pinjong(p):
    pol, (r, rg, ra, rf) = tannprofil(p)
    rb = p["boring"] / 2
    z4, z5 = p["b"], p["b"] + p["l_nav"]
    m = mf.Manifold.extrude(mf.CrossSection([pol]), p["b"]) + syl(p["d_nav"] / 2, z4, z5)
    m = m - syl(rb, -1, z5 + 1)
    m = m - syl(rb + 0.4, -0.01, 0.4, rb) - syl(rb, z5 - 0.4, z5 + 0.01, rb + 0.4)  # innføringsfaser
    return m - klemnav_kutt(p, rb, z4, z5), (r, rg, ra, rf)


def tegn(p, m, radier, sti):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    r, rg, ra, rf = radier
    z4, z5 = p["b"], p["b"] + p["l_nav"]
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, z, tittel in ((a, p["b"] / 2, "Tannhjulssnitt (midt i tennene)"), (b, (z4 + z5) / 2, "Navsnitt med skruene")):
        for poly in m.slice(z).to_polygons():
            ax.add_patch(Polygon(np.asarray(poly), fc="#d4952f", ec="k", lw=0.6))
        ax.set_title(tittel); ax.set_xlim(-11.5, 11.5); ax.set_ylim(-11.5, 11.5); ax.set_aspect("equal")
    for rr, navn in ((r, "delesirkel"), (rg, "grunnsirkel")):
        a.add_patch(plt.Circle((0, 0), rr, fc="none", ec="#2a5db5", lw=0.6, ls="--"))
    a.text(0, -10.8, f"Ø{2 * ra:.1f} topp / Ø{2 * r:.1f} deling / Ø{2 * rf:.1f} rot", ha="center", fontsize=8)
    snitt = m.rotate([90, 0, 0]).slice(0.0)
    for poly in snitt.to_polygons():
        c.add_patch(Polygon(np.asarray(poly) * [1, -1], fc="#d4952f", ec="k", lw=0.6))
    c.set_title("Snitt langs aksen (y = 0)"); c.set_xlabel("x [mm]"); c.set_ylabel("z [mm] (0 = mot svinghjulet)")
    c.set_xlim(-11.5, 11.5); c.set_ylim(-1, z5 + 1); c.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(sti, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slakk", type=float, default=P["slakk"])
    ap.add_argument("--boring", type=float, default=P["boring"])
    a = ap.parse_args()
    p = dict(P, slakk=a.slakk, boring=a.boring)
    UT.mkdir(parents=True, exist_ok=True)
    m, radier = pinjong(p)
    assert m.status() == mf.Error.NoError, "ugyldig geometri"
    navn = f"pinjong_z{p['z']}_m{p['modul']:g}"
    skriv_3mf(m.rotate([180, 0, 0]).translate([0, 0, p["b"] + p["l_nav"]]), UT / f"{navn}.3mf", f"Pinjong z{p['z']} m{p['modul']:g} med klemnav")
    tegn(p, m, radier, UT / "pinjong_forhandsvisning.png")
    r, rg, ra, rf = radier
    print(f"{navn}: Ø topp {2 * ra:.2f}, Ø deling {2 * r:.2f}, Ø rot {2 * rf:.2f}, "
          f"vegg rot-boring {rf - p['boring'] / 2:.2f}, lengde {p['b'] + p['l_nav']:.1f}, "
          f"y i bilen 5 til {5 + p['b'] + p['l_nav']:.1f}, masse {m.volume() * 1.24e-3:.2f} g")
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
