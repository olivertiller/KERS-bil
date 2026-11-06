"""Figur av den generative rammen til rapporten, laget fra STL-filen som ble printet.

    python3 Code/rammefigur.py

Skriver to riss til Rapport/Images/Ramme/: ramme_skraatt.png (skrått forfra) og ramme_side.png
(fra siden). Hvitt rundt modellen klippes bort.
Skriver også ut ytre mål, volum og masse i PLA, som er tallene i rapporten.
"""

import struct
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image, ImageChops

ROT = Path(__file__).resolve().parent.parent
STL = ROT / "Bygg" / "Ramme v1.stl"
UT = ROT / "Rapport" / "Images" / "Ramme"
RHO_PLA = 1.24e-3  # g/mm^3


def les_stl(fil):
    data = fil.read_bytes()
    n = struct.unpack("<I", data[80:84])[0]
    dt = np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")])
    return np.frombuffer(data[84 : 84 + 50 * n], dtype=dt)["v"].reshape(-1, 3, 3).astype(float)


def skygge(tri, lys):
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    k = 0.35 + 0.65 * np.clip(n @ lys, 0, 1)
    farge = np.array([0.62, 0.66, 0.70])
    return np.column_stack([np.outer(k, farge), np.ones(len(k))])


def tegn(ax, tri, elev, azim):
    # Lyset følger kameraet, så flatene mot betrakteren blir lyse.
    e, a = np.radians(elev), np.radians(azim)
    kamera = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    lys = kamera + np.array([0.0, 0.0, 0.6])
    lys /= np.linalg.norm(lys)
    ax.add_collection3d(Poly3DCollection(tri, facecolors=skygge(tri, lys), edgecolor="none", antialiased=False))
    mn, mx = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
    pad = 0.03 * (mx - mn).max()
    ax.set_xlim(mn[0] - pad, mx[0] + pad); ax.set_ylim(mn[1] - pad, mx[1] + pad); ax.set_zlim(mn[2] - pad, mx[2] + pad)
    ax.set_box_aspect(mx - mn + 2 * pad, zoom=0.95)
    ax.set_proj_type("ortho")
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()


def main():
    tri = les_stl(STL)
    mn, mx = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
    volum = np.sum(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))) / 6
    print(f"Ytre mål: {' x '.join(f'{d:.1f}' for d in mx - mn)} mm")
    print(f"Volum: {volum / 1000:.1f} cm^3, massiv PLA: {volum * RHO_PLA:.0f} g")

    UT.mkdir(parents=True, exist_ok=True)
    for navn, elev, azim in [("ramme_skraatt", 28, -58), ("ramme_side", 0, -90)]:
        fig = plt.figure(figsize=(6, 6), dpi=300)
        tegn(fig.add_subplot(projection="3d"), tri, elev, azim)
        fil = UT / f"{navn}.png"
        fig.savefig(fil, facecolor="white")
        plt.close(fig)
        klipp(fil)
        print(f"Skrev {fil.relative_to(ROT)}")


def klipp(fil, marg=20):
    bilde = Image.open(fil).convert("RGB")
    boks = ImageChops.difference(bilde, Image.new("RGB", bilde.size, "white")).getbbox()
    v, o, h, n = boks
    bilde.crop((max(v - marg, 0), max(o - marg, 0), h + marg, n + marg)).save(fil)

if __name__ == "__main__":
    main()
