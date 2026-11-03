"""Kontroll og forhåndsvisning av oppsettet i KERS_GD_oppsett.py, uten Fusion.

  * Sjekker at ingen bevar-kropp overlapper en hindring (da avviser Fusion studien).
  * Regner volum og PLA-masse for bevar-kroppene og lokkene.
  * Tegner sideriss, planriss og frontriss til Kers-bil/Generativ_ramme/gd_oppsett.png.

    python3 Code/Fusion/gd_forhandsvisning.py
"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

HER = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("gd", HER / "KERS_GD_oppsett" / "KERS_GD_oppsett.py")
gd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gd)
UT = HER.parent.parent / "Kers-bil" / "Generativ_ramme"
RHO_PLA = 1.24e-3  # g/mm^3


def inni(prim, X, Y, Z):
    t = prim[0]
    if t == "boks":
        _, x0, x1, y0, y1, z0, z1 = prim
        return (X > x0) & (X < x1) & (Y > y0) & (Y < y1) & (Z > z0) & (Z < z1)
    if t == "syl_y":
        _, x, z, r, y0, y1 = prim
        return ((X - x) ** 2 + (Z - z) ** 2 < r * r) & (Y > y0) & (Y < y1)
    _, x, y, r, z0, z1 = prim
    return ((X - x) ** 2 + (Y - y) ** 2 < r * r) & (Z > z0) & (Z < z1)


def inni_kropp(k, X, Y, Z):
    m = np.zeros(X.shape, bool)
    for tegn, p in k["ops"]:
        m = (m | inni(p, X, Y, Z)) if tegn == "+" else (m & ~inni(p, X, Y, Z))
    return m


def bbox(prim):
    t = prim[0]
    if t == "boks":
        return prim[1:]
    if t == "syl_y":
        _, x, z, r, y0, y1 = prim
        return (x - r, x + r, y0, y1, z - r, z + r)
    _, x, y, r, z0, z1 = prim
    return (x - r, x + r, y - r, y + r, z0, z1)


def punkter(k, h=0.25):
    b = [bbox(p) for s, p in k["ops"] if s == "+"]
    lo = np.min([[q[0], q[2], q[4]] for q in b], axis=0)
    hi = np.max([[q[1], q[3], q[5]] for q in b], axis=0)
    ax = [np.arange(lo[i] + h / 2 + 1e-3, hi[i], h) for i in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    m = inni_kropp(k, X, Y, Z)
    return X[m], Y[m], Z[m], h**3


def kontroll(kropper):
    bevar = [k for k in kropper if k["gruppe"] == "Bevar"]
    hindr = [k for k in kropper if k["gruppe"] in ("Hindring", "Grenser")]
    feil = 0
    print("Bevar-kropper (volum, masse i PLA):")
    tot = 0.0
    for k in bevar:
        X, Y, Z, dv = punkter(k)
        v = len(X) * dv
        tot += v
        print(f"  {k['navn']:<24} {v:7.0f} mm³  {v * RHO_PLA:5.2f} g")
        for o in hindr:
            n = int(inni_kropp(o, X, Y, Z).sum())
            if n:
                feil += 1
                print(f"    OVERLAPP med {o['navn']}: {n * dv:.1f} mm³")
    print(f"  {'Sum':<24} {tot:7.0f} mm³  {tot * RHO_PLA:5.2f} g")
    for k in kropper:
        if k["gruppe"] == "Lagerlokk":
            X, Y, Z, dv = punkter(k)
            print(f"  {k['navn']:<24} {len(X) * dv:7.0f} mm³  {len(X) * dv * RHO_PLA:5.2f} g  (printes separat)")
    print("Ingen overlapp mellom bevar og hindring." if feil == 0 else f"{feil} overlapp må fjernes!")
    return feil


def tegn(kropper):
    farger = {"Bevar": ("#2e8b57", 0.75), "Hindring": ("#c0392b", 0.12), "Lagerlokk": ("#e0a030", 0.8)}
    fig, akser = plt.subplots(3, 1, figsize=(9, 15))
    # (tittel, akse for u, akse for v) der 1 = x, 3 = y, 5 = z (indeks i boks-tuppelen)
    riss = [("Sideriss (x–z), fra +y", 1, 5, (-12, 152), (-5, 82)),
            ("Planriss (x–y), ovenfra", 1, 3, (-12, 152), (-52, 52)),
            ("Frontriss (y–z), bakfra", 3, 5, (-52, 52), (-5, 82))]
    for ax, (tittel, a, b, xl, yl) in zip(akser, riss):
        ax.set_title(tittel)
        for gruppe in ("Hindring", "Lagerlokk", "Bevar"):
            fc, al = farger[gruppe]
            for k in kropper:
                if k["gruppe"] != gruppe or k["navn"].startswith("Grense"):
                    continue
                for tegn_, p in k["ops"]:
                    if tegn_ != "+":
                        continue
                    t = p[0]
                    if t == "boks":
                        bb = p[1:]
                        u0, u1 = (bb[0], bb[1]) if a == 1 else (bb[2], bb[3])
                        v0, v1 = (bb[2], bb[3]) if b == 3 else (bb[4], bb[5])
                        ax.add_patch(Rectangle((u0, v0), u1 - u0, v1 - v0, fc=fc, ec=fc, alpha=al, lw=.6))
                    elif (t == "syl_y" and (a, b) == (1, 5)) or (t == "syl_z" and (a, b) == (1, 3)):
                        c = (p[1], p[2])
                        ax.add_patch(Circle(c, p[3], fc=fc, ec=fc, alpha=al, lw=.6))
                    else:
                        bb = bbox(p)
                        u0, u1 = (bb[0], bb[1]) if a == 1 else (bb[2], bb[3])
                        v0, v1 = (bb[2], bb[3]) if b == 3 else (bb[4], bb[5])
                        ax.add_patch(Rectangle((u0, v0), u1 - u0, v1 - v0, fc=fc, ec=fc, alpha=al, lw=.6))
        if a == 1 and b == 5:
            ax.add_patch(Rectangle((0, 0), 140, 70, fill=False, ls="--", ec="k", lw=.8))
        elif b == 3:
            ax.add_patch(Rectangle((0, -42.5), 140, 85, fill=False, ls="--", ec="k", lw=.8))
        else:
            ax.add_patch(Rectangle((-42.5, 0), 85, 70, fill=False, ls="--", ec="k", lw=.8))
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
    akser[0].text(2, 76, "grønt = bevar   rødt = hindring   oransje = lagerlokk   stiplet = 140 × 85 × 70",
                  fontsize=8)
    fig.tight_layout()
    UT.mkdir(parents=True, exist_ok=True)
    fig.savefig(UT / "gd_oppsett.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    k = gd.geometri()
    print(f"{len(k)} kropper: " + ", ".join(f"{g} {sum(1 for x in k if x['gruppe'] == g)}"
                                           for g in ("Bevar", "Hindring", "Grenser", "Lagerlokk")))
    kontroll(k)
    tegn(k)
    print(f"Bilde: {UT / 'gd_oppsett.png'}")
