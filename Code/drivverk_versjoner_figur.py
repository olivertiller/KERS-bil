"""Setter sammen sidebildene av de fire drivverkversjonene til én figur på en rekke, med markeringer.

Bildene lages av Code/Fusion/KERS_drivverk_versjoner (drivverk_v<n>_side.png). Målestokken finnes fra
underlaget i bildet, som går fra x = -15 til 160 mm og har overflaten på z = 0. Sett fra tannhjulssiden
går x mot venstre i bildet.

    python3 Code/drivverk_versjoner_figur.py
gir Rapport/Images/Drivverk/drivverk_versjoner.pdf
"""

import argparse
import importlib.util
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from PIL import Image

REPO = Path(__file__).resolve().parent.parent
BILDER = REPO / "Rapport" / "Images" / "Drivverk"
_spec = importlib.util.spec_from_file_location(
    "dv", REPO / "Code" / "Fusion" / "KERS_drivverk_versjoner" / "KERS_drivverk_versjoner.py")
dv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dv)

UNDERLAG = (-15.0, 160.0)
UTSNITT = (-3.0, 146.0, -3.0, 71.0)  # x0, x1, z0, z1 [mm]
ROD, GRA = "#d62728", "#555555"
LYSGRA = "#c8c8c8"  # markering av uendret klaring, synlig mot det mørke svinghjulet
TITLER = {1: "V1: konsept 3", 2: "V2: aksler i lik høyde", 3: "V3: Ø6-aksler", 4: "V4: $a$ = 36 mm"}
FARGE_UNDERLAG = np.array([201, 200, 192])


def last(v):
    return np.asarray(Image.open(BILDER / f"drivverk_v{v}_side.png").convert("RGB")).astype(int)


def fjern_rutenett(im):
    """Rutenettet (nøytralt grått, 205-250) og origo-aksene (rød/grønn/blå, 2 kanaler like) blir hvite."""
    im = im.copy()
    mx, mn = im.max(axis=2), im.min(axis=2)
    noytral = (mx - mn < 4) & (mn >= 205)
    srt = np.sort(im, axis=2)
    akse = (srt[..., 2] >= 220) & (np.abs(srt[..., 1] - srt[..., 0]) < 6) & (srt[..., 0] >= 110) & (srt[..., 0] <= 160)
    im[noytral | akse] = 255
    return im


def malestokk(im, underlag):
    """Finn px per mm og pikselposisjonen til x = 0 og z = 0 fra underlaget."""
    treff = np.all(np.abs(im - FARGE_UNDERLAG) < 6, axis=2)
    rader = np.where(treff.sum(axis=1) > im.shape[1] * 0.5)[0]
    if len(rader) == 0:
        raise SystemExit("Fant ikke underlaget i bildet")
    rad = rader[len(rader) // 2]
    kol = np.where(treff[rad])[0]
    venstre, hoyre = kol.min() - 1, kol.max() + 1  # svart kantlinje
    sk = (hoyre - venstre) / (underlag[1] - underlag[0])
    topp = rader.min() - 1  # kantlinjen over er z = 0
    return lambda x, z: (hoyre - (x - underlag[0]) * sk, topp - z * sk), sk


def klaring_punkt(p):
    """Midt i gapet mellom svinghjulskanten og forakselen."""
    dx, dz = p["x_for"] - p["x_sv"], p["z_hjul"] - p["z_sv"]
    d = math.hypot(dx, dz)
    t = (dv.R_SV + p["klaring_aksel"] / 2) / d
    return p["x_sv"] + dx * t, p["z_sv"] + dz * t


def tekst(ax, xy, s, farge, **kw):
    ax.annotate(s, xy=xy, fontsize=7, color=farge, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), **kw)


def tegn(ax, v, underlag):
    im = last(v)
    T, sk = malestokk(im, underlag)
    im = fjern_rutenett(im)
    ox, oz = T(0, 0)
    yy, xx = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    im[((xx - ox) ** 2 + (yy - oz) ** 2 < (3 * sk) ** 2) & (yy < oz - 1)] = 255
    p = dv.mal(v)
    x0, x1, z0, z1 = UTSNITT
    (u0, w1), (u1, w0) = T(x1, z1), T(x0, z0)
    ax.imshow(im[int(w1):int(w0), int(u0):int(u1)].astype(np.uint8),
              extent=(u0, u1, w0, w1))
    ax.set_xlim(u0, u1)
    ax.set_ylim(w0, w1)
    ax.axis("off")
    ax.set_title(TITLER[v], fontsize=8, pad=2)

    # Klaringen mellom svinghjul og foraksel, i alle versjonene
    kx, kz = klaring_punkt(p)
    endret = v in (3, 4)
    farge = ROD if endret else GRA
    ax.add_patch(Circle(T(kx, kz), 6 * sk, fill=False, ec=ROD if endret else LYSGRA, lw=1.4 if endret else 1.1))
    mm = f"{p['klaring_aksel']:.1f}".replace(".", ",").replace(",0", "")
    tekst(ax, T(kx - 15, kz + 15), f"{mm} mm", farge)

    if v == 1:
        # Ulik høyde på akslene
        for x, z, dx in ((p["x_bak"], p["z_hjul"], -8), (p["x_sv"], p["z_sv"], 6)):
            ax.plot(*zip(T(x - 11, z), T(x + 11, z)), color=GRA, lw=0.9, ls="--")
            tekst(ax, T(x + dx, z - 11), f"{z:.0f} mm", GRA)
    if v == 2:
        # Alle akslene i lik høyde, større hjul
        ax.plot(*zip(T(p["x_bak"] - 8, 34), T(p["x_for"] + 8, 34)), color=ROD, lw=1.2, ls="--")
        ax.annotate("", xy=T(p["x_for"], 0.3), xytext=T(p["x_for"], 2 * p["r_hjul"] - 0.3),
                    arrowprops=dict(arrowstyle="<->", color=ROD, lw=1.0))
        tekst(ax, T(p["x_for"], 9), "Ø68", ROD)
    if v == 3:
        ax.annotate("Ø6-aksel", xy=T(p["x_for"], p["z_hjul"] - 3.5), xytext=T(p["x_for"], 10),
                    fontsize=7, color=ROD, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                    arrowprops=dict(arrowstyle="->", color=ROD, lw=1.0))
    if v == 4:
        ax.add_patch(Circle(T(p["x_bak"], p["z_hjul"]), (dv.M * 75 / 2 + 2) * sk, fill=False, ec=ROD, lw=1.2))
        tekst(ax, T(p["x_bak"], p["z_hjul"] + 14), "73 tenner", ROD)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammelt-underlag", action="store_true",
                    help="bilder fra første kjøring, der underlaget gikk til lengden + 15")
    ap.add_argument("--ut", default=str(BILDER / "drivverk_versjoner.pdf"))
    a = ap.parse_args()
    fig, axs = plt.subplots(1, 4, figsize=(6.3, 1.05), dpi=300)
    for ax, v in zip(axs, (1, 2, 3, 4)):
        und = (UNDERLAG[0], dv.mal(v)["lengde"] + 15) if a.gammelt_underlag else UNDERLAG
        tegn(ax, v, und)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.86, bottom=0.0, wspace=0.03)
    fig.savefig(a.ut, dpi=300)
    print("skrev", a.ut)


if __name__ == "__main__":
    main()
