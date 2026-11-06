"""Maskintegninger (A4 liggende, førstevinkelprojeksjon) av snortrommel TR-1 og drivhjul DH-2.

    python3 Code/maskintegninger.py

Geometrien hentes fra de samme parameterne som 3MF-filene (snortrommel.py og hjul.py), og snitt og
riss lages fra modellen med manifold3d, så tegningen viser nøyaktig det som printes. Krever
manifold3d og matplotlib. Skriver Images/Tegninger/TR-1.pdf og Images/Tegninger/DH-2.pdf.

Arket er i mm (x mot høyre, y opp). Hvert riss har sin egen avbildning fra modellen til arket,
valgt etter ISO 128 metode E (førstevinkel): riss sett fra høyre står til venstre, sett ovenfra
står under.
"""

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

import hjul as H
import snortrommel as T

UT = Path(__file__).resolve().parent.parent / "Images" / "Tegninger"
TYKK, TYNN = 0.5, 0.25  # linjebredder i mm
PT = 72 / 25.4  # pt per mm
H_MAL, H_TITTEL = 2.5, 5.0  # teksthøyder i mm
PIL = 2.5  # pillengde i mm
STREKPUNKT = (0, (14, 3, 2, 3))
STIPLET = (0, (4, 2))
plt.rcParams["hatch.linewidth"] = 0.25 * PT


def tall(x):
    """Desimalkomma, uten unødvendige desimaler."""
    return f"{x:.2f}".rstrip("0").rstrip(".").replace(".", ",")


class Ark:
    """Et A4-ark i mm med ramme og tittelfelt."""

    def __init__(self):
        self.fig = plt.figure(figsize=(297 / 25.4, 210 / 25.4))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, 297)
        self.ax.set_ylim(0, 210)
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        self.ax.add_patch(Rectangle((20, 10), 267, 190, fill=False, lw=0.7 * PT))

    # --- grunnelementer ---
    def linje(self, pts, lw=TYKK, ls="-", z=3):
        pts = np.asarray(pts, float)
        self.ax.plot(pts[:, 0], pts[:, 1], color="k", lw=lw * PT, ls=ls, solid_capstyle="butt",
                     dash_capstyle="butt", zorder=z)

    def tekst(self, x, y, s, ha="center", va="baseline", rot=0, h=H_MAL, **kw):
        self.ax.text(x, y, s, ha=ha, va=va, rotation=rot, fontsize=h / 0.72 * PT, family="DejaVu Sans",
                     rotation_mode="anchor", zorder=6, **kw)

    def senterlinje(self, p0, p1):
        self.linje([p0, p1], lw=TYNN, ls=STREKPUNKT, z=2)

    def kryss(self, c, r):
        """Senterlinjekryss for en sirkel med radius r (arkets mm)."""
        self.senterlinje(np.asarray(c) - [r + 3, 0], np.asarray(c) + [r + 3, 0])
        self.senterlinje(np.asarray(c) - [0, r + 3], np.asarray(c) + [0, r + 3])

    def sirkel(self, c, r, lw=TYKK, ls="-"):
        th = np.linspace(0, 2 * np.pi, 181)
        self.linje(np.c_[c[0] + r * np.cos(th), c[1] + r * np.sin(th)], lw=lw, ls=ls)

    def flate(self, polys, skravur=None, klipp=None):
        """Polygoner med kantlinjer. skravur: None, '45' eller '-45' (snittflate). Hull blir åpne.
        klipp = (senter, radius) begrenser alt til en sirkel (detaljriss)."""
        verts, koder = [], []
        for p in polys:
            verts += list(p) + [p[0]]
            koder += [MPath.MOVETO] + [MPath.LINETO] * (len(p) - 1) + [MPath.CLOSEPOLY]
        sti = MPath(verts, koder)
        artister = [self.ax.add_patch(PathPatch(sti, fc="white", ec="none", zorder=3))]
        if skravur:
            artister.append(self.ax.add_patch(PathPatch(sti, fc="none", ec="k", lw=0,
                                              hatch="////" if skravur == "45" else "\\\\\\\\", zorder=3)))
        for p in polys:
            q = np.vstack([p, p[:1]])
            artister += self.ax.plot(q[:, 0], q[:, 1], color="k", lw=TYKK * PT, zorder=4)
        if klipp:
            sirkel = plt.Circle(klipp[0], klipp[1], transform=self.ax.transData, fc="none", ec="none")
            self.ax.add_patch(sirkel)
            for a in artister:
                a.set_clip_path(sirkel)

    def detaljmerke(self, c, r, bokstav, vinkel=45):
        """Tynn sirkel rundt et område som vises forstørret, med bokstav."""
        self.sirkel(c, r, lw=TYNN)
        v = math.radians(vinkel)
        p = np.asarray(c) + (r + 3.5) * np.array([math.cos(v), math.sin(v)])
        self.tekst(p[0], p[1] - 1.8, bokstav, h=H_TITTEL)

    def iso(self, m, x0, y0, w, h, elev=22, azim=-58, rot=None):
        """Skyggelagt isometrisk riss av modellen i et felt (mm), uten målestokk."""
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        mesh = m.to_mesh() if rot is None else m.rotate(rot).to_mesh()
        v = np.asarray(mesh.vert_properties)[:, :3]
        tri = v[np.asarray(mesh.tri_verts)]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
        lys = np.array([-0.4, -0.6, 0.7])
        lys /= np.linalg.norm(lys)
        g = 0.45 + 0.5 * np.clip(n @ lys, 0, 1)
        ax3 = self.fig.add_axes([x0 / 297, y0 / 210, w / 297, h / 210], projection="3d")
        ax3.add_collection3d(Poly3DCollection(tri, facecolors=np.c_[g, g, g * 1.02 - 0.02, np.ones_like(g)],
                                              edgecolors="none"))
        lo, hi = v.min(0), v.max(0)
        c, r = (lo + hi) / 2, (hi - lo).max() / 2
        ax3.set_xlim(c[0] - r, c[0] + r); ax3.set_ylim(c[1] - r, c[1] + r); ax3.set_zlim(c[2] - r, c[2] + r)
        ax3.set_box_aspect((1, 1, 1))
        ax3.view_init(elev=elev, azim=azim)
        ax3.set_proj_type("ortho")
        ax3.axis("off")
        ax3.patch.set_alpha(0)
        ax3.set_facecolor((1, 1, 1, 0))

    def pil(self, spiss, retning):
        d = np.asarray(retning, float)
        d = d / np.linalg.norm(d)
        n = np.array([-d[1], d[0]])
        b = np.asarray(spiss, float) - d * PIL
        self.ax.add_patch(Polygon([spiss, b + n * PIL / 6, b - n * PIL / 6], fc="k", ec="none", zorder=6))

    # --- mål ---
    def _maalelinje(self, a, b, tekst, rot, side, ute=1):
        """Målelinje a-b med piler og tekst. side = +1 tekst over/til venstre."""
        a, b = np.asarray(a, float), np.asarray(b, float)
        d = (b - a) / np.linalg.norm(b - a)
        L = np.linalg.norm(b - a)
        n = np.array([-d[1], d[0]])
        bredde = 0.62 * H_MAL * len(tekst) + 1
        if L >= 2 * PIL + 1:
            self.linje([a, b], lw=TYNN)
            self.pil(a, -d)
            self.pil(b, d)
            m = (a + b) / 2
            if L < bredde:  # teksten får ikke plass over målelinjen, legg den utenfor
                m = (b if ute > 0 else a) + ute * d * (bredde / 2 + 2)
                self.linje([b, b + d * (bredde + 3)] if ute > 0 else [a, a - d * (bredde + 3)], lw=TYNN)
        else:  # pilene utenfor
            if ute > 0:
                self.linje([a - d * (PIL + 3), b + d * (PIL + bredde + 3)], lw=TYNN)
                m = b + d * (PIL + 1.5 + bredde / 2)
            else:
                self.linje([a - d * (PIL + bredde + 3), b + d * (PIL + 3)], lw=TYNN)
                m = a - d * (PIL + 1.5 + bredde / 2)
            self.pil(a, d)
            self.pil(b, -d)
        self.tekst(*(m + n * 0.9), tekst, rot=rot, va="baseline")

    def maal_h(self, p1, p2, y, tekst, ute=1):
        """Vannrett mål mellom x-ene til p1 og p2, målelinje i høyden y. Hjelpelinjer fra punktene."""
        for p in (p1, p2):
            s = np.sign(y - p[1]) or 1
            if abs(y - p[1]) > 1.2:
                self.linje([[p[0], p[1] + s * 1.0], [p[0], y + s * 2.0]], lw=TYNN)
        x1, x2 = sorted((p1[0], p2[0]))
        self._maalelinje([x1, y], [x2, y], tekst, 0, 1, ute)

    def maal_v(self, p1, p2, x, tekst, ute=1):
        """Loddrett mål mellom y-ene til p1 og p2, målelinje ved x. Teksten leses fra høyre."""
        for p in (p1, p2):
            s = np.sign(x - p[0]) or 1
            if abs(x - p[0]) > 1.2:
                self.linje([[p[0] + s * 1.0, p[1]], [x + s * 2.0, p[1]]], lw=TYNN)
        y1, y2 = sorted((p1[1], p2[1]))
        self._maalelinje([x, y1], [x, y2], tekst, 90, 1, ute)

    def henvisning(self, fra, til, tekst, ha="left", prikk=False):
        """Henvisningslinje fra delen (pil, eller prikk på en flate) til tekst med vannrett skulder."""
        fra, til = np.asarray(fra, float), np.asarray(til, float)
        dx = 5 if ha == "left" else -5
        self.linje([fra, til, til + [dx, 0]], lw=TYNN)
        if prikk:
            self.ax.add_patch(plt.Circle(fra, 0.5, color="k", zorder=6))
        else:
            self.pil(fra, fra - til)
        self.tekst(til[0] + dx + (0.8 if ha == "left" else -0.8), til[1] - 0.9, tekst, ha=ha)

    def snittlinje(self, p, retning, bokstav, ende):
        """Ende av en snittlinje: tykk strek ved p, pil i synsretningen og bokstaven."""
        p, d, e = (np.asarray(v, float) for v in (p, retning, ende))
        self.linje([p, p + e * 5], lw=0.7)
        q = p + e * 5
        self.linje([q, q + d * 7], lw=TYNN)
        self.pil(q + d * 7, d)
        self.tekst(*(q + d * 5 + e * 3.2 - [0, H_TITTEL / 2.6]), bokstav, h=H_TITTEL)

    # --- tittelfelt ---
    def tittelfelt(self, tittel, nr, materiale, masse, malestokk, merknad):
        x0, y0, w = 107, 10, 180
        for y in (7, 14, 21):
            self.linje([[x0, y0 + y], [x0 + w, y0 + y]], lw=0.35)
        self.linje([[x0, y0 + 32], [x0 + w, y0 + 32]], lw=0.7)
        self.linje([[x0, y0], [x0, y0 + 32]], lw=0.7)
        for x in (60, 120):
            self.linje([[x0 + x, y0], [x0 + x, y0 + 21]], lw=0.35)
        felter = [
            (0, 0, "Tegnet", "Oliver Kvaløy Tiller"), (60, 0, "Dato", "2026-11-06"), (120, 0, "Ark", "1 av 1"),
            (0, 7, "Materiale", materiale), (60, 7, "Masse", masse), (120, 7, "Målestokk", malestokk),
            (0, 14, "Generelle toleranser", "ISO 2768-m"), (60, 14, "Mål", "mm"), (120, 14, "Tegningsnr.", nr),
        ]
        for x, y, led, verdi in felter:
            self.tekst(x0 + x + 1.5, y0 + y + 4.8, led, ha="left", h=1.8, color="0.35")
            self.tekst(x0 + x + 1.5, y0 + y + 1.4, verdi, ha="left", h=2.5)
        self.tekst(x0 + 1.5, y0 + 28.6, "NTNU · MAST1101 · Svinghjulsdrevet lekebil", ha="left", h=1.8, color="0.35")
        self.tekst(x0 + 1.5, y0 + 22.5, tittel, ha="left", h=H_TITTEL, weight="bold")
        self.projeksjonssymbol(x0 + w - 16, y0 + 26.5)
        for i, s in enumerate(merknad):
            self.tekst(25, 14 + 4.5 * (len(merknad) - 1 - i), s, ha="left", h=2.5)

    def projeksjonssymbol(self, x, y):
        """Førstevinkelprojeksjon (ISO 128 metode E): kjegle til venstre, endebilde til høyre."""
        h = 3.0
        self.linje([[x - 10, y - h / 2], [x - 10, y + h / 2], [x - 4, y + h], [x - 4, y - h], [x - 10, y - h / 2]],
                   lw=0.35)
        for r in (h, h / 2):
            self.sirkel((x + 4, y), r, lw=0.35)
        self.senterlinje([x - 12, y], [x + 8.5, y])

    def lagre(self, sti):
        self.fig.savefig(sti)
        plt.close(self.fig)


def polygoner(cs):
    return [np.asarray(p, float) for p in cs.to_polygons()]


def avbild(polys, o, k, f):
    """f: (n, 2) modellkoordinater i snittplanet -> (n, 2) risskoordinater (før skala)."""
    return [np.asarray(o) + k * f(p) for p in polys]


# ----------------------------------------------------------------------------------------------
# TR-1
# ----------------------------------------------------------------------------------------------

def tegn_tr1():
    p, d_k = T.P, 12.0
    m = T.trommel(p, d_k)
    rk, rf, (z1, z2, z3, z4, z5) = T.lengder(p, d_k)
    rb, rn, ru = p["boring"] / 2, p["d_nav"] / 2, p["d_utsparing"] / 2
    zm = (z4 + z5) / 2
    zs, ze = z1 + p["krok_klaring"], z1 + p["krok_klaring"] + p["krok_bredde"]
    s_, L_, g_, h_ = p["krok_stilk"] / 2, p["krok_finger"], p["krok_gap"], p["krok_hoyde"]
    e, xs = p["e_skrue"], p["x_sete"]
    k = 2.0  # målestokk 2:1
    ende = lambda q: np.c_[-q[:, 1], q[:, 0]]  # snitt vinkelrett på aksen, sett fra navenden: (x, y) -> (-y, x)
    ark = Ark()

    # Snitt A-A (hovedriss): planet y = 0. Arket: u = z (navet til høyre), v = x (kroken ned).
    oA = np.array([122.0, 136.0])
    A = lambda z, x: oA + k * np.array([z, x])
    s_A = m.rotate([90, 0, 0]).slice(0.0)  # (x, y, z) -> (x, -z, y): planet y = 0 med koordinater (x, -z)
    ark.flate(avbild(polygoner(s_A), oA, k, lambda q: np.c_[-q[:, 1], q[:, 0]]), skravur="45")
    ark.senterlinje(A(-3, 0), A(z5 + 3, 0))
    ark.tekst(*A(z5 / 2, rn + 17), "A–A", h=H_TITTEL)
    y1, y2 = A(0, rn + 4)[1], A(0, rn + 8.5)[1]
    ark.maal_h(A(0, rf), A(z1, rf), y1, tall(z1))
    ark.maal_h(A(z1, rf), A(z2, rk), y1, tall(z2 - z1))
    ark.maal_h(A(z2, rk), A(z4, rf), y1, tall(z4 - z2))
    ark.maal_h(A(z4, rn), A(z5, rn), y1, tall(z5 - z4))
    ark.maal_h(A(0, rf), A(z5, rn), y2, tall(z5))
    yk = A(0, -(rn + 6))[1]
    ark.maal_h(A(z1, -rf), A(zs, -rk), yk, tall(zs - z1), ute=-1)
    ark.maal_h(A(zs, -rk), A(ze, -(rk + h_)), yk, tall(ze - zs))
    ark.maal_v(A(0, -ru), A(0, ru), A(-4, 0)[0], "Ø" + tall(2 * ru))
    ark.maal_v(A(0, -rf), A(0, rf), A(-8.5, 0)[0], "Ø" + tall(2 * rf))
    ark.maal_v(A(z5, -rb), A(z5, rb), A(z5 + 4, 0)[0], "Ø" + tall(2 * rb))
    ark.maal_v(A(z5, -rn), A(z5, rn), A(z5 + 8.5, 0)[0], "Ø" + tall(2 * rn))
    ark.henvisning(A((z2 + z3) / 2 + 0.3, (rk + rf) / 2 + 0.1), A(z2 + 5, rn + 2.2), "45°")
    ark.henvisning(A(1.6, ru - 1.2), A(3.5, 1.2), "45°")
    # Snittlinjer B-B (gjennom kroken) og C-C (gjennom skruene), sett mot svinghjulet (-z)
    zb = 7.5
    for z, bokstav in ((zb, "B"), (zm, "C")):
        r0 = rn + 9.5
        ark.snittlinje(A(z, r0), [-1, 0], bokstav, [0, 1])
        ark.snittlinje(A(z, -r0), [-1, 0], bokstav, [0, -1])
        ark.linje([A(z, r0), A(z, -r0)], lw=TYNN, ls=STREKPUNKT, z=5)

    # Snitt C-C gjennom skrueaksene, til venstre. Arket: u = -y, v = x.
    oC = np.array([55.0, 136.0])
    C = lambda x, y: oC + k * np.array([-y, x])
    ark.flate(avbild(polygoner(m.slice(zm)), oC, k, ende), skravur="45")
    ark.kryss(oC, k * rn)
    ark.tekst(*(oC + [0, k * rn + 25]), "C–C", h=H_TITTEL)
    ark.maal_h(C(rn - 1, e), C(rn - 1, -e), C(rn + 9, 0)[1], tall(2 * e))
    ark.maal_v(C(0, -(e + p["mutter"] / 2)), C(-xs, -(e + p["mutter"] / 2)), C(0, -(rn + 3))[0], tall(xs))
    ark.maal_h(C(-rn + 0.5, e - p["d_hode"] / 2), C(-rn + 0.5, e + p["d_hode"] / 2), C(-(rn + 8), 0)[1],
               "Ø" + tall(p["d_hode"]))
    ark.maal_h(C(-rn + 0.5, -e - p["mutter"] / 2), C(-rn + 0.5, -e + p["mutter"] / 2), C(-(rn + 8), 0)[1],
               tall(p["mutter"]))
    ark.henvisning(C(-2.2, e + p["m3_hull"] / 2), np.array([44.0, 112.0]), f"Ø{tall(p['m3_hull'])}", ha="right")
    ark.henvisning(C(0, rn - 1.2) + [0, 0.6], np.array([44.0, 166.0]), f"Slisse {tall(p['slisse'])}", ha="right")
    ark.tekst(*(C(-(rn + 8), -e) + [0, -6.5]), "mutter", h=2.0, color="0.25")
    ark.tekst(*(C(-(rn + 8), e) + [0, -6.5]), "skruehode", h=2.0, color="0.25")
    # Snittlinje A-A: planet y = 0 er den loddrette senterlinjen; A-A ses fra +y, dvs. mot +u
    ark.snittlinje(oC + [0, k * rn + 2], [1, 0], "A", [0, 1])
    ark.snittlinje(oC + [0, -(k * rn + 2)], [1, 0], "A", [0, -1])

    # Snitt B-B gjennom kroken, til høyre. Arket: u = -y, v = x.
    oB = np.array([222.0, 136.0])
    ark.flate(avbild(polygoner(m.slice(zb)), oB, k, ende), skravur="45")
    ark.kryss(oB, k * rf)
    ark.sirkel(oB, k * rf, lw=TYNN, ls=STIPLET)  # flensen bak, skjult kant
    ark.tekst(*(oB + [0, k * rf + 9]), "B–B", h=H_TITTEL)
    ark.maal_v(oB + k * np.array([0, -rk]), oB + k * np.array([0, rk]), oB[0] + k * (rf + 4), "Ø" + tall(2 * rk))
    K = lambda t, r, o=oB, kk=k: o + kk * np.array([t, -(rk + r)])  # krokens (t, r) -> arket
    cE = K(-s_ - L_ / 2 + 0.5, h_ / 2)
    ark.detaljmerke(cE, k * 3.6, "E", vinkel=200)
    th = np.linspace(math.radians(25), math.radians(75), 30)  # opptrekk: mot klokka sett fra navenden
    rr = k * rf + 5
    bue = np.c_[oB[0] + rr * np.cos(th), oB[1] + rr * np.sin(th)]
    ark.linje(bue, lw=TYNN)
    ark.pil(bue[-1], bue[-1] - bue[-3])
    ark.tekst(oB[0] + rr * 0.93 + 1, oB[1] + rr * 0.33, "Opptrekk", ha="left")

    # Detalj E: kroken, 5:1
    kd = 5.0
    tE, rE = -s_ - L_ / 2 + 0.5, h_ / 2
    oE = np.array([258.0, 90.0]) - kd * np.array([tE, -(rk + rE)])
    E = lambda t, r: K(t, r, oE, kd)
    rkl = 4.6
    ark.flate(avbild(polygoner(m.slice(zb)), oE, kd, ende), skravur="45", klipp=(E(tE, rE), kd * rkl))
    ark.sirkel(E(tE, rE), kd * rkl, lw=TYNN)
    ark.tekst(*(E(tE, rE) + [0, kd * rkl + 4]), "E (5:1)", h=H_TITTEL)
    yb = E(0, h_ + 1.2)[1]
    ark.maal_h(E(-s_ - L_, h_), E(-s_, g_), yb, tall(L_))
    ark.maal_h(E(-s_, g_), E(s_, h_), yb, tall(2 * s_))
    ark.maal_v(E(s_, 0), E(s_, h_), E(s_ + 1.8, 0)[0], tall(h_))
    ark.maal_v(E(-s_ - L_, 0.04), E(-s_ - L_, g_), E(-s_ - L_ - 1.4, 0)[0], tall(g_))
    ark.henvisning(E(-s_ - L_ * 0.6, g_), E(-s_ - L_ - 1.2, -0.6) + [-4, 0], "Løkka", ha="right", prikk=True)

    # Isometrisk riss nede til venstre
    ark.iso(m, 26, 18, 66, 56, elev=24, azim=-50, rot=[0, 0, 0])
    ark.tekst(59, 16, "Isometrisk, uten målestokk", h=2.0, color="0.35")

    masse = m.volume() * T.RHO_PLA
    merknad = ["Utsparing for sveisekransen: Ø14,5, 1 mm rett, deretter 45° inn til boringen.",
               "Printes stående med indre flens mot platen, uten support. Klemnavet strammes med",
               "2 × M3 × 12 og mutter. Løkka på snora legges rundt kroken under fingeren og glir av",
               "av seg selv når snora er trukket helt ut."]
    for i, s in enumerate(merknad):
        ark.tekst(110, 60 - 4.5 * i, s, ha="left", h=2.5)
    ark.tittelfelt("Snortrommel TR-1", "KERS-TR-1", "PLA", f"{tall(round(masse, 1))} g", "2:1 (E 5:1)", [])
    ark.lagre(UT / "TR-1.pdf")
    return masse


# ----------------------------------------------------------------------------------------------
# DH-2: felg i PLA (pos. 1) og dekk i TPU (pos. 2)
# ----------------------------------------------------------------------------------------------

def tegn_dh2():
    p = H.P
    felg, dekk = H.hjul(p, True), H.dekk(p)
    b = p["bredde"]
    r_sete = p["d_ytre"] / 2 - p["t_dekk"]
    rs = r_sete + p["skulder"]
    r_inn = r_sete - p["t_felg"]
    kant = (b - p["b_dekk"]) / 2
    rn, rkr, rr_ = p["d_nav"] / 2, p["d_krage"] / 2, p["d_rund"] / 2
    zt = b + p["l_krage"]
    ts = p["t_dekk"]
    bd, f = p["b_dekk"] - 0.2, 0.4
    zd = kant + 0.1  # dekket montert midt mellom skuldrene
    ro = r_sete + ts  # montert dekk, ytre radius
    k = 1.0  # målestokk 1:1
    ark = Ark()
    sideriss = lambda q: np.c_[-q[:, 1], q[:, 0]]  # (x, z) -> arket (-z, x)

    def snitt_y0(mm):
        """Planet y = 0 som (x, z)-polygoner."""
        return [np.c_[q[:, 0], -q[:, 1]] for q in polygoner(mm.rotate([90, 0, 0]).slice(0.0))]

    def dekkprofil(r0, z0):
        return np.array([(z0, r0), (z0, r0 + ts - f), (z0 + f, r0 + ts), (z0 + bd - f, r0 + ts),
                         (z0 + bd, r0 + ts - f), (z0 + bd, r0)])

    def tegn_dekk(o, kk, r0, z0, sgn=1, klipp=None):
        pr = dekkprofil(r0, z0)
        ark.flate([o + kk * np.c_[-pr[:, 0], sgn * pr[:, 1]]], skravur="-45", klipp=klipp)

    # Sett fra utsiden (ser mot +z): arket (u, v) = (y, x), så en eike står loddrett.
    oF = np.array([68.0, 142.0])
    ark.flate(avbild(polygoner(felg.project()), oF, k, lambda q: q[:, ::-1]))
    ark.sirkel(oF, k * ro)  # dekket utenfor skulderen
    ark.kryss(oF, k * ro)
    a = math.radians(60)
    ark.henvisning(oF + k * 21 * np.array([math.sin(a), math.cos(a)]), oF + [32, 28],
                   f"{p['n_eiker']} eiker × {tall(p['b_eike'])}")
    ark.detaljmerke(oF, k * rn + 2.5, "C", vinkel=-45)
    ark.snittlinje(oF + [0, k * ro + 3], [1, 0], "A", [0, 1])
    ark.snittlinje(oF + [0, -(k * ro + 3)], [1, 0], "A", [0, -1])

    # Snitt A-A til høyre: ses i +y-retning, arket (u, v) = (-z, x). Utsiden (z = 0) til høyre.
    oA = np.array([150.0, 142.0])
    A = lambda z, x: oA + k * np.array([-z, x])
    ark.flate(avbild(snitt_y0(felg), oA, k, sideriss), skravur="45")
    for sgn in (1, -1):
        tegn_dekk(oA, k, r_sete, zd, sgn)
    ark.senterlinje(A(zt + 3, 0), A(-3, 0))
    ark.tekst(*A(b / 2, ro + 15), "A–A", h=H_TITTEL)
    ark.maal_h(A(0, rs), A(b, rs), A(0, ro + 9)[1], tall(b))
    ark.maal_h(A(b, -rn), A(zt, -rkr), A(0, -(rn + 6))[1], tall(p["l_krage"]))
    ark.maal_h(A(0, -(rn + 8)), A(p["t_skive"], -(rn + 8)), A(0, -(rn + 12))[1], tall(p["t_skive"]))
    for i, (r, z, navn) in enumerate([(rn, 0, f"Ø{tall(2 * rn)}"), (r_inn, 0, f"Ø{tall(2 * r_inn)}"),
                                      (ro, zd + bd / 2, f"Ø{tall(2 * ro)}")]):
        ark.maal_v(A(z, -r), A(z, r), A(-5 - 6.5 * i, 0)[0], navn)
    ark.maal_v(A(zt, -rkr), A(zt, rkr), A(zt + 5, 0)[0], "Ø" + tall(2 * rkr))
    ark.henvisning(A(zt - 0.6, rr_ - 0.3), A(zt + 7, rn + 6), f"Ø{tall(2 * rr_)}, dybde {tall(p['l_rund'])}",
                   ha="right")
    ark.henvisning(A(b * 0.45, ro - 0.5), A(b + 8, ro + 6), "Pos. 2 dekk", ha="right", prikk=True)
    ark.henvisning(A(b * 0.75, r_inn + 1), A(b + 8, r_inn - 9), "Pos. 1 felg", ha="right", prikk=True)
    ark.detaljmerke(A(b / 2, r_sete + 0.3), 7.0, "D", vinkel=-30)

    # Detalj D: felgen med dekket, 5:1
    kd = 5.0
    cD_mod = (b / 2, r_sete - 0.4)
    oD = np.array([75.0, 58.5]) - kd * np.array([-cD_mod[0], cD_mod[1]])
    D = lambda z, x: oD + kd * np.array([-z, x])
    klipp = (D(*cD_mod), kd * 4.4)
    ark.flate(avbild(snitt_y0(felg), oD, kd, sideriss), skravur="45", klipp=klipp)
    tegn_dekk(oD, kd, r_sete, zd, klipp=klipp)
    ark.sirkel(*klipp, lw=TYNN)
    ark.tekst(*(klipp[0] + [0, kd * 4.4 + 4]), "D (5:1)", h=H_TITTEL)
    ark.maal_h(D(b, rs), D(b - kant, rs), D(0, ro + 1.4)[1], tall(kant))
    ark.maal_h(D(b - kant, r_sete), D(kant, r_sete), D(0, ro + 1.4)[1], tall(b - 2 * kant))
    ark.maal_v(D(0, r_inn), D(0, r_sete), D(-1.6, 0)[0], tall(p["t_felg"]))
    ark.maal_v(D(0, r_sete), D(0, rs), D(-3.2, 0)[0], tall(p["skulder"]))
    ark.henvisning(D(kant + 0.6, r_sete), D(-6.5, r_sete - 3.0), f"Felgsete Ø{tall(2 * r_sete)}", prikk=True)

    # Detalj C: firkanthullet, 5:1, sett fra utsiden
    kc = 5.0
    oC = np.array([240.0, 156.0])
    for q in polygoner(felg.slice(4.0)):
        if np.max(np.abs(q)) < p["firkant"]:
            ark.linje(np.vstack([oC + kc * q[:, ::-1], oC + kc * q[:1, ::-1]]))
    ark.sirkel(oC, kc * rn)
    ark.sirkel(oC, kc * rkr, lw=TYNN, ls=STIPLET)  # kragen på innsiden, skjult
    ark.kryss(oC, kc * p["firkant"] / 2 + 1)
    s2 = p["firkant"] / 2
    ark.maal_h(oC + kc * np.array([-s2, -s2]), oC + kc * np.array([s2, -s2]), oC[1] - kc * s2 - 7, tall(p["firkant"]))
    ark.maal_v(oC + kc * np.array([s2, -s2]), oC + kc * np.array([s2, s2]), oC[0] + kc * s2 + 7, tall(p["firkant"]))
    c = s2 - 0.45 * 0.35
    ark.henvisning(oC + kc * (np.array([-c, c]) + 0.45 * np.array([-0.71, 0.71])), oC + [-18, kc * s2 + 10],
                   "4 × R0,45", ha="right")
    ark.tekst(*(oC + [0, kc * rn + 5]), "C (5:1)", h=H_TITTEL)

    # Pos. 2 dekk før montering, snitt 5:1 (bare en side)
    kt = 5.0
    oT = np.array([268.0, 92.0])
    P2 = lambda z, r: oT + kt * np.array([-z, r])
    tegn_dekk(oT, kt, 0.0, 0.0)
    ark.maal_h(P2(0, ts), P2(bd, ts), P2(0, ts + 2.4)[1], tall(bd))
    ark.maal_v(P2(0, 0), P2(0, ts), P2(-1.6, 0)[0], tall(ts))
    ark.henvisning(P2(bd - f / 2, ts - f / 2), P2(bd + 1.4, ts + 1.2), f"{tall(f)} × 45°", ha="right")
    ark.tekst(*(P2(bd / 2, ts + 3.6)), "Pos. 2 før montering (5:1)", h=3.5)
    rdi = r_sete - p["strekk"] / 2
    ark.tekst(*(P2(bd / 2, -2.6)), f"Innvendig Ø{tall(2 * rdi)}, strekkes på Ø{tall(2 * r_sete)}", h=2.5)

    # Isometrisk riss
    iso_m = felg + dekk.translate([0, 0, zd]).scale([r_sete / (r_sete - p["strekk"] / 2)] * 2 + [1])
    ark.iso(iso_m, 150, 44, 78, 64, elev=38, azim=-35)
    ark.tekst(189, 46, "Isometrisk, uten målestokk", h=2.0, color="0.35")

    m_felg, m_dekk = felg.volume() * H.RHO_PLA, dekk.volume() * H.RHO_TPU
    merknad = ["Pos. 1 felg i PLA, printes med utsiden ned.", "Pos. 2 dekk i TPU, strekkes over felgsetet.",
               "Akselenden freses til 4 × 4 mm.", "Kragen ligger an mot lagerets innerring."]
    for i, s in enumerate(merknad):
        ark.tekst(25, 32 - 4.5 * i, s, ha="left", h=2.5)
    ark.tittelfelt("Drivhjul DH-2", "KERS-DH-2", "PLA (1), TPU (2)",
                   f"{tall(round(m_felg, 1))} g + {tall(round(m_dekk, 1))} g", "1:1 (C, D 5:1)", [])
    ark.lagre(UT / "DH-2.pdf")
    return m_felg, m_dekk


def main():
    UT.mkdir(parents=True, exist_ok=True)
    m = tegn_tr1()
    mf_, md = tegn_dh2()
    print(f"TR-1: {m:.1f} g, DH-2: felg {mf_:.1f} g, dekk {md:.1f} g. Skrevet til {UT}")


if __name__ == "__main__":
    main()
