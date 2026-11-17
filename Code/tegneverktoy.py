"""Verktøy for maskintegninger på A4 i målestokk: riss og snitt med synlige og skjulte kanter
beregnet fra 3D-modellen (manifold3d), mål, henvisninger, posisjonsnumre, stykkliste og tittelfelt.

Arket er i mm med origo nede til venstre, og PDF-en lagres i full størrelse (297 x 210 eller
210 x 297), slik at et mål på arket delt på målestokken gir delens mål.

Riss: en modell vises med en basis (R, U): R peker mot høyre på arket, U opp, og betrakteren ser
langs D = U x R (inn i arket). Kanter som tegnes er skarpe kanter (vinkel mellom flatene over
grensen) og konturkanter (der flaten snur fra betrakteren). Hver kant deles i punkter, og et punkt
er synlig hvis en stråle fra punktet mot betrakteren ikke treffer modellen. Snitt fjerner delen
foran snittplanet og skraverer snittflaten.
"""

import math

import manifold3d as mf
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

TYKK, TYNN = 0.5, 0.25  # linjebredder i mm (ISO 128: 0,5 / 0,25)
PT = 72 / 25.4  # pt per mm
H_MAL, H_TITTEL = 2.5, 5.0
PIL = 2.5
STREKPUNKT = (0, (12, 3, 1.5, 3))
STIPLET = (0, (3, 1.5))
plt.rcParams["hatch.linewidth"] = 0.25 * PT
plt.rcParams["pdf.fonttype"] = 42
MASKE = dict(fc="white", ec="none", pad=0.4)  # hvit bakgrunn bak tekst, så linjer aldri krysser teksten


def tall(x):
    """Desimalkomma, uten unødvendige desimaler."""
    return f"{x:.2f}".rstrip("0").rstrip(".").replace(".", ",")


# ------------------------------------------------------------------------------------------
# Modeller
# ------------------------------------------------------------------------------------------

def fra_trekanter(v, t):
    """Manifold fra hjørner og trekanter. Like hjørner slås sammen."""
    v = np.asarray(v, np.float32)
    nokkel = np.round(v.astype(np.float64), 4)
    _, inv = np.unique(nokkel, axis=0, return_inverse=True)
    inv = inv.ravel()
    unike = np.zeros((inv.max() + 1, 3), np.float32)
    unike[inv] = v
    t = inv[np.asarray(t)]
    t = t[(t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])]
    m = mf.Manifold(mf.Mesh(vert_properties=unike, tri_verts=t.astype(np.uint32)))
    assert m.status() == mf.Error.NoError, m.status()
    return m


def les_stl(sti):
    import struct

    b = open(sti, "rb").read()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:84 + n * 50], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
    v = a["v"].reshape(-1, 3)
    return fra_trekanter(v, np.arange(len(v)).reshape(-1, 3))


# ------------------------------------------------------------------------------------------
# Skjult-linje-beregning
# ------------------------------------------------------------------------------------------

def basis(R, U):
    R, U = np.asarray(R, float), np.asarray(U, float)
    W = np.cross(R, U)  # mot betrakteren
    return np.array([R, U, W])


def _kjed(segmenter, tol=1e-3):
    """Setter sammen segmenter (a, b) i 2D til polylinjer gjennom punkter med to naboer."""
    if not segmenter:
        return []
    S = np.asarray(segmenter)
    P = np.round(S.reshape(-1, 2) / tol).astype(np.int64)
    _, idx = np.unique(P, axis=0, return_inverse=True)
    idx = idx.ravel().reshape(-1, 2)
    naboer = {}
    for k, (a, b) in enumerate(idx):
        if a == b:
            continue
        naboer.setdefault(a, []).append((k, b))
        naboer.setdefault(b, []).append((k, a))
    koord = {}
    for k, (a, b) in enumerate(idx):
        koord[a], koord[b] = S[k, 0], S[k, 1]
    brukt = np.zeros(len(idx), bool)
    linjer = []

    def gaa(start, k0):
        sti, cur, k = [start], start, k0
        while True:
            brukt[k] = True
            a, b = idx[k]
            nxt = b if a == cur else a
            sti.append(nxt)
            cur = nxt
            videre = [(kk, o) for kk, o in naboer.get(cur, []) if not brukt[kk]]
            if len(naboer.get(cur, [])) != 2 or not videre:
                return sti
            k = videre[0][0]

    for node, liste in naboer.items():  # start i endepunkter og forgreninger
        if len(liste) != 2:
            for k, _ in liste:
                if not brukt[k]:
                    linjer.append(gaa(node, k))
    for k in range(len(idx)):  # lukkede sløyfer
        if not brukt[k] and idx[k, 0] != idx[k, 1]:
            linjer.append(gaa(idx[k, 0], k))
    return [np.array([koord[i] for i in l]) for l in linjer]


def til_vis(m, B, snitt=None):
    """Modellen i rissets koordinater (x = R, y = U, z mot betrakteren), uten det som ligger foran
    snittplanet. Returnerer (modell, snittflater)."""
    mv = m.transform(np.c_[B, np.zeros(3)])
    flater = []
    if snitt is not None:
        flater = [np.asarray(p, float) for p in mv.slice(snitt - 1e-4).to_polygons()]
        mv = mv.trim_by_plane([0, 0, -1], -snitt)
    return mv, flater


def kanter(m, B, vinkel=20.0, snitt=None, skjulte=False, steg=0.15, tol=1e-6, vis=None, okkluder=None):
    """Synlige (og evt. skjulte) kanter for modellen m sett med basis B.
    snitt: dybde (langs W) for snittplanet; alt foran (nærmere betrakteren) fjernes.
    vis: ferdig (modell, snittflater) fra til_vis. okkluder: modell i rissets koordinater som kan
    skjule kantene (i en sammenstilling alle delene), ellers delen selv.
    Returnerer (synlige, skjulte, snittflater) i rissets koordinater (x = R, y = U)."""
    mv, flater = vis if vis is not None else til_vis(m, B, snitt)
    occ = okkluder if okkluder is not None else mv
    mesh = mv.to_mesh()
    V = np.asarray(mesh.vert_properties, float)[:, :3]
    T = np.asarray(mesh.tri_verts, np.int64)
    if len(T) == 0:
        return [], [], flater
    n = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])
    A = np.linalg.norm(n, axis=1)
    n = n / (A[:, None] + 1e-30)
    # kantliste med flatene på hver side
    e = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    f = np.tile(np.arange(len(T)), 3)
    es = np.sort(e, axis=1)
    orden = np.lexsort((es[:, 1], es[:, 0]))
    es, f = es[orden], f[orden]
    lik = np.all(es[1:] == es[:-1], axis=1)
    i1 = np.nonzero(lik)[0]
    ea, f1, f2 = es[i1], f[i1], f[i1 + 1]
    cosv = np.sum(n[f1] * n[f2], axis=1)
    skarp = cosv < math.cos(math.radians(vinkel))
    foran1, foran2 = n[f1, 2] > 1e-6, n[f2, 2] > 1e-6
    kontur = foran1 != foran2
    sett = skarp | kontur
    if not skjulte:
        sett &= foran1 | foran2 | kontur
    ea = ea[sett]
    p0, p1 = V[ea[:, 0]], V[ea[:, 1]]
    L = np.linalg.norm((p1 - p0)[:, :2], axis=1)
    keep = L > 1e-5
    p0, p1, L = p0[keep], p1[keep], L[keep]
    zmax = max(V[:, 2].max(), occ.bounding_box()[5]) + 1.0
    eps = 1e-3
    synlige, skjulte_seg = [], []
    for a, b, l in zip(p0, p1, L):
        k = max(2, int(math.ceil(l / steg)) + 1)
        t = np.linspace(0, 1, k)
        mid = (t[:-1] + t[1:]) / 2
        punkter = a + (b - a) * mid[:, None]
        sy = np.empty(len(mid), bool)
        for i, q in enumerate(punkter):
            treff = occ.ray_cast([q[0], q[1], q[2] + eps], [q[0], q[1], zmax])
            sy[i] = not any(h.distance > tol for h in treff)
        # del kanten i løp med lik synlighet
        grenser = np.nonzero(np.diff(sy.astype(int)))[0]
        start = 0
        for g in list(grenser) + [len(sy) - 1]:
            ta, tb = t[start], t[g + 1]
            seg = (a[:2] + (b - a)[:2] * ta, a[:2] + (b - a)[:2] * tb)
            (synlige if sy[start] else skjulte_seg).append(seg)
            start = g + 1
    return _kjed(synlige), (_kjed(skjulte_seg) if skjulte else []), flater


# ------------------------------------------------------------------------------------------
# Arket
# ------------------------------------------------------------------------------------------

class Ark:
    """Et A4-ark i mm med ramme. liggende: 297 x 210, ellers 210 x 297."""

    def __init__(self, liggende=True):
        self.B, self.H = (297, 210) if liggende else (210, 297)
        self.fig = plt.figure(figsize=(self.B / 25.4, self.H / 25.4))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, self.B)
        self.ax.set_ylim(0, self.H)
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        self.ax.add_patch(Rectangle((20, 10), self.B - 30, self.H - 20, fill=False, lw=0.7 * PT))
        # sentreringsmerker (ISO 5457)
        for x, y, dx, dy in ((self.B / 2 + 5, 10, 0, -5), (self.B / 2 + 5, self.H - 10, 0, 5),
                             (20, self.H / 2, -5, 0), (self.B - 10, self.H / 2, 5, 0)):
            self.linje([[x, y], [x + dx, y + dy]], lw=0.7)

    # --- grunnelementer ---
    def linje(self, pts, lw=TYKK, ls="-", z=3, farge="k"):
        pts = np.asarray(pts, float)
        self.ax.plot(pts[:, 0], pts[:, 1], color=farge, lw=lw * PT, ls=ls, solid_capstyle="butt",
                     dash_capstyle="butt", solid_joinstyle="round", zorder=z)

    def tekst(self, x, y, s, ha="center", va="baseline", rot=0, h=H_MAL, maske=False, **kw):
        if maske:
            kw["bbox"] = MASKE
        self.ax.text(x, y, s, ha=ha, va=va, rotation=rot, fontsize=h / 0.72 * PT, family="DejaVu Sans",
                     rotation_mode="anchor", **{"zorder": 7, **kw})

    def senterlinje(self, p0, p1):
        self.linje([p0, p1], lw=TYNN, ls=STREKPUNKT, z=2)

    def kryss(self, c, r, ut=3):
        """Senterlinjekryss for en sirkel med radius r (arkets mm)."""
        c = np.asarray(c, float)
        self.senterlinje(c - [r + ut, 0], c + [r + ut, 0])
        self.senterlinje(c - [0, r + ut], c + [0, r + ut])

    def sirkel(self, c, r, lw=TYKK, ls="-", bue=(0, 360), z=3):
        th = np.radians(np.linspace(bue[0], bue[1], 181))
        self.linje(np.c_[c[0] + r * np.cos(th), c[1] + r * np.sin(th)], lw=lw, ls=ls, z=z)

    def flate(self, polys, skravur="45", tetthet=4, z=3):
        """Skravert flate (snittflate). skravur: '45', '-45', '45t' (tett) eller None (hvit)."""
        verts, koder = [], []
        for p in polys:
            verts += list(p) + [p[0]]
            koder += [MPath.MOVETO] + [MPath.LINETO] * (len(p) - 1) + [MPath.CLOSEPOLY]
        sti = MPath(verts, koder)
        self.ax.add_patch(PathPatch(sti, fc="white", ec="none", zorder=z))
        if skravur:
            monster = {"45": "/" * tetthet, "-45": "\\" * tetthet}[skravur]
            self.ax.add_patch(PathPatch(sti, fc="none", ec="k", lw=0, hatch=monster, zorder=z))

    def pil(self, spiss, retning):
        d = np.asarray(retning, float)
        d = d / np.linalg.norm(d)
        n = np.array([-d[1], d[0]])
        b = np.asarray(spiss, float) - d * PIL
        self.ax.add_patch(Polygon([spiss, b + n * PIL / 6, b - n * PIL / 6], fc="k", ec="none", zorder=6))

    def prikk(self, p, r=0.5):
        self.ax.add_patch(plt.Circle(p, r, color="k", zorder=6))

    # --- mål ---
    def _maalelinje(self, a, b, tekst, rot, ute=1, tekst_pos=None):
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
            if L < bredde:  # teksten får ikke plass, legg den utenfor
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
        if tekst_pos is not None:
            m = a + d * tekst_pos
        self.tekst(*(m + n * 0.9), tekst, rot=rot, va="baseline", maske=True)

    def maal_h(self, p1, p2, y, tekst, ute=1, tekst_pos=None):
        """Vannrett mål mellom x-ene til p1 og p2, målelinje i høyden y."""
        for p in (p1, p2):
            s = np.sign(y - p[1]) or 1
            if abs(y - p[1]) > 1.2:
                self.linje([[p[0], p[1] + s * 1.0], [p[0], y + s * 2.0]], lw=TYNN)
        x1, x2 = sorted((p1[0], p2[0]))
        self._maalelinje([x1, y], [x2, y], tekst, 0, ute, tekst_pos)

    def maal_v(self, p1, p2, x, tekst, ute=1, tekst_pos=None):
        """Loddrett mål mellom y-ene til p1 og p2, målelinje ved x. Teksten leses fra høyre."""
        for p in (p1, p2):
            s = np.sign(x - p[0]) or 1
            if abs(x - p[0]) > 1.2:
                self.linje([[p[0] + s * 1.0, p[1]], [x + s * 2.0, p[1]]], lw=TYNN)
        y1, y2 = sorted((p1[1], p2[1]))
        self._maalelinje([x, y1], [x, y2], tekst, 90, ute, tekst_pos)

    def maal_dia(self, c, r, vinkel, tekst, ute=0.0):
        """Diametermål tvers over en sirkel (arkets mm) i retning vinkel (grader)."""
        c = np.asarray(c, float)
        d = np.array([math.cos(math.radians(vinkel)), math.sin(math.radians(vinkel))])
        a, b = c - d * r, c + d * r
        self.linje([a, b + d * ute], lw=TYNN)
        self.pil(a, -d)
        self.pil(b, d)
        rot = vinkel if -90 < vinkel <= 90 else vinkel - 180
        m = c + d * (r * 0.45 if ute == 0 else r + ute * 0.6)
        n = np.array([-d[1], d[0]]) if -90 < vinkel <= 90 else np.array([d[1], -d[0]])
        self.tekst(*(m + n * 0.9), tekst, rot=rot, maske=True)

    def maal_radius(self, c, r, vinkel, tekst, ut=8):
        """Radiusmål: linje fra sentrum gjennom buen med pil på buen, tekst utenfor."""
        c = np.asarray(c, float)
        d = np.array([math.cos(math.radians(vinkel)), math.sin(math.radians(vinkel))])
        b = c + d * r
        self.linje([c, b + d * ut], lw=TYNN)
        self.pil(b, d)
        e = b + d * ut
        self.linje([e, e + [4 if d[0] >= 0 else -4, 0]], lw=TYNN)
        self.tekst(e[0] + (5 if d[0] >= 0 else -5), e[1] - 0.9, tekst, ha="left" if d[0] >= 0 else "right")

    def vinkelmaal(self, c, r, a0, a1, tekst, hjelpelinjer=()):
        """Vinkelmål: bue med sentrum c og radius r fra a0 til a1 (grader), piler i begge ender."""
        c = np.asarray(c, float)
        th = np.radians(np.linspace(a0, a1, 60))
        bue = np.c_[c[0] + r * np.cos(th), c[1] + r * np.sin(th)]
        self.linje(bue, lw=TYNN)
        t0 = np.array([-math.sin(th[0]), math.cos(th[0])])
        t1 = np.array([-math.sin(th[-1]), math.cos(th[-1])])
        self.pil(bue[0], -t0)
        self.pil(bue[-1], t1)
        for p0, p1 in hjelpelinjer:
            self.linje([p0, p1], lw=TYNN)
        am = math.radians((a0 + a1) / 2)
        p = c + (r + 2.5) * np.array([math.cos(am), math.sin(am)])
        self.tekst(p[0], p[1] - 1.0, tekst, maske=True)

    def henvisning(self, fra, til, tekst, ha="left", prikk=False, pil=True):
        """Henvisningslinje fra delen (pil på en kant, prikk på en flate) til tekst med skulder."""
        fra, til = np.asarray(fra, float), np.asarray(til, float)
        dx = 4 if ha == "left" else -4
        self.linje([fra, til, til + [dx, 0]], lw=TYNN)
        if prikk:
            self.prikk(fra)
        elif pil:
            self.pil(fra, fra - til)
        self.tekst(til[0] + dx + (0.8 if ha == "left" else -0.8), til[1] - 0.9, tekst, ha=ha)

    def ballong(self, fra, sentrum, nr, r=4.0):
        """Posisjonsnummer i sirkel med henvisningslinje som ender i en prikk på delen."""
        s = np.asarray(sentrum, float)
        if fra is not None:  # uten henvisningslinje: gruppert ved siden av en annen ballong
            fra = np.asarray(fra, float)
            d = fra - s
            d = d / np.linalg.norm(d)
            self.linje([s + d * r, fra], lw=TYNN, z=6)
            self.prikk(fra)
        self.ax.add_patch(plt.Circle(s, r, fc="white", ec="k", lw=TYNN * PT, zorder=6))
        self.tekst(s[0], s[1] - 1.75, str(nr), h=3.5, zorder=8)

    def snittlinje(self, p, retning, bokstav, ende):
        """Ende av en snittlinje: tykk strek ved p, pil i synsretningen og bokstaven."""
        p, d, e = (np.asarray(v, float) for v in (p, retning, ende))
        self.linje([p, p + e * 5], lw=0.7)
        q = p + e * 5
        self.linje([q, q + d * 7], lw=TYNN)
        self.pil(q + d * 7, d)
        self.tekst(*(q + d * 5 + e * 3.2 - [0, H_TITTEL / 2.6]), bokstav, h=H_TITTEL)

    def snittplan(self, a, b, retning, bokstav):
        """Snittlinje fra a til b (arkets mm): tynn strekpunkt, tykke ender, piler i synsretningen."""
        a, b = np.asarray(a, float), np.asarray(b, float)
        e = (b - a) / np.linalg.norm(b - a)
        self.linje([a, b], lw=TYNN, ls=STREKPUNKT, z=5)
        self.snittlinje(b, retning, bokstav, e)
        self.snittlinje(a, retning, bokstav, -e)

    def detaljmerke(self, c, r, bokstav, vinkel=45):
        self.sirkel(c, r, lw=TYNN)
        v = math.radians(vinkel)
        p = np.asarray(c) + (r + 3.5) * np.array([math.cos(v), math.sin(v)])
        self.tekst(p[0], p[1] - 1.8, bokstav, h=H_TITTEL)

    def overskrift(self, x, y, s, h=H_TITTEL):
        self.tekst(x, y, s, h=h)

    # --- tittelfelt og stykkliste ---
    def tittelfelt(self, tittel, nr, materiale, masse, malestokk, type_="Arbeidstegning", dato="2026-11-17",
                   rev="B"):
        x0, y0, w = self.B - 10 - 180, 10, 180
        self.tf = (x0, y0 + 32)
        for y in (7, 14, 21):
            self.linje([[x0, y0 + y], [x0 + w, y0 + y]], lw=0.35)
        self.linje([[x0, y0 + 32], [x0 + w, y0 + 32]], lw=0.7)
        self.linje([[x0, y0], [x0, y0 + 32]], lw=0.7)
        for x in (60, 120):
            self.linje([[x0 + x, y0], [x0 + x, y0 + 21]], lw=0.35)
        self.linje([[x0 + 150, y0], [x0 + 150, y0 + 7]], lw=0.35)
        felter = [
            (0, 0, "Tegnet", "Oliver Kvaløy Tiller"), (60, 0, "Dato", dato), (120, 0, "Ark", "1 av 1"),
            (150, 0, "Rev.", rev),
            (0, 7, "Materiale", materiale), (60, 7, "Masse", masse), (120, 7, "Målestokk", malestokk),
            (0, 14, "Generelle toleranser", "ISO 2768-m"), (60, 14, "Dokumenttype", type_),
            (120, 14, "Tegningsnr.", nr),
        ]
        for x, y, led, verdi in felter:
            self.tekst(x0 + x + 1.5, y0 + y + 4.8, led, ha="left", h=1.8, color="0.35")
            self.tekst(x0 + x + 1.5, y0 + y + 1.4, verdi, ha="left", h=2.5)
        self.tekst(x0 + 1.5, y0 + 28.6, "NTNU · MAST1101 · Svinghjulsdrevet lekebil · mål i mm", ha="left", h=1.8,
                   color="0.35")
        self.tekst(x0 + 1.5, y0 + 22.5, tittel, ha="left", h=H_TITTEL, weight="bold")
        self.projeksjonssymbol(x0 + w - 16, y0 + 26.5)

    def stykkliste(self, rader, bredder=(10, 12, 58, 50, 50), h=6.0):
        """Stykkliste rett over tittelfeltet, nederste rad er overskriften (ISO 7573).
        rader: (pos, antall, betegnelse, materiale/standard, tegningsnr./merknad)."""
        x0, y0 = self.tf
        w = sum(bredder)
        x0 = self.B - 10 - w
        hode = ("Pos.", "Ant.", "Betegnelse", "Materiale / standard", "Tegningsnr. / merknad")
        alle = [hode] + list(rader)
        for i, rad in enumerate(alle):
            y = y0 + i * h
            self.linje([[x0, y + h], [x0 + w, y + h]], lw=0.35 if i else 0.7)
            x = x0
            for j, (b, s) in enumerate(zip(bredder, rad)):
                self.tekst(x + (b / 2 if j < 2 else 1.5), y + 1.9, str(s), ha="center" if j < 2 else "left",
                           h=2.2 if i else 2.0, color="k" if i else "0.35")
                x += b
        yt = y0 + len(alle) * h
        x = x0
        for b in (0,) + tuple(bredder):
            x += b
            self.linje([[x, y0], [x, yt]], lw=0.35 if 0 < x - x0 < w else 0.7)
        self.linje([[x0, yt], [x0 + w, yt]], lw=0.7)
        return yt

    def tabell(self, x0, y0, rader, bredder, h=5.0, tittel=None):
        """Enkel tabell med øvre venstre hjørne i (x0, y0)."""
        w = sum(bredder)
        n = len(rader) + (1 if tittel else 0)
        for i in range(n + 1):
            self.linje([[x0, y0 - i * h], [x0 + w, y0 - i * h]], lw=0.7 if i in (0, n) else 0.35)
        rad0 = 0
        if tittel:
            self.tekst(x0 + 1.5, y0 - h + 1.6, tittel, ha="left", h=2.5, weight="bold")
            rad0 = 1
        x = x0
        for j, b in enumerate((0,) + tuple(bredder)):
            x += b
            self.linje([[x, y0 - rad0 * h], [x, y0 - n * h]], lw=0.7 if j in (0, len(bredder)) else 0.35)
        self.linje([[x0, y0], [x0, y0 - n * h]], lw=0.7)
        self.linje([[x0 + w, y0], [x0 + w, y0 - n * h]], lw=0.7)
        for i, rad in enumerate(rader):
            x = x0
            for b, s in zip(bredder, rad):
                self.tekst(x + 1.5, y0 - (i + rad0 + 1) * h + 1.6, str(s), ha="left", h=2.3)
                x += b

    def projeksjonssymbol(self, x, y):
        """Førstevinkelprojeksjon (ISO 128 metode E): kjegle til venstre, endebilde til høyre."""
        h = 3.0
        self.linje([[x - 10, y - h / 2], [x - 10, y + h / 2], [x - 4, y + h], [x - 4, y - h], [x - 10, y - h / 2]],
                   lw=0.35)
        for r in (h, h / 2):
            self.sirkel((x + 4, y), r, lw=0.35)
        self.senterlinje([x - 12, y], [x + 8.5, y])

    def merknader(self, x, y, linjer, h=2.5):
        for i, s in enumerate(linjer):
            self.tekst(x, y - 4.2 * i, s, ha="left", h=h)

    def lagre(self, sti):
        self.fig.savefig(sti)
        plt.close(self.fig)


# ------------------------------------------------------------------------------------------
# Riss
# ------------------------------------------------------------------------------------------

class Riss:
    """Et riss eller snitt av modellen m på arket.

    R, U: modellretninger mot høyre og opp på arket. o: arkpunktet der modellpunktet ref havner.
    snitt: et modellpunkt i snittplanet (planet står vinkelrett på synsretningen).
    """

    def __init__(self, ark, m, R, U, o, k, ref=(0, 0, 0), snitt=None, skjulte=False, skravur="45",
                 vinkel=20.0, steg=None, tegn=True, tol=1e-6):
        self.ark, self.k = ark, k
        self.B = basis(R, U)
        self.o = np.asarray(o, float)
        self.ref = np.asarray(ref, float)
        self.m = m
        d_snitt = None if snitt is None else float(self.B[2] @ np.asarray(snitt, float))
        steg = steg or max(0.04, 0.6 / k)
        self.syn, self.skj, self.flater = kanter(m, self.B, vinkel, d_snitt, skjulte, steg, tol)
        self.skravur = skravur
        if tegn:
            self.tegn()

    def a(self, uv):
        """Rissets koordinater (modellmm langs R, U) -> arket."""
        uv = np.asarray(uv, float)
        r0 = self.B[:2] @ self.ref
        return self.o + self.k * (uv - r0)

    def p(self, x, y=None, z=None):
        """Modellpunkt -> arket."""
        q = np.asarray(x if y is None else (x, y, z), float)
        return self.o + self.k * (self.B[:2] @ (q - self.ref))

    def tegn(self, flater=True, klipp=None):
        """klipp = (arkpunkt, radius): bare det som ligger innenfor sirkelen tegnes (detaljriss)."""
        ark = self.ark
        n0 = len(ark.ax.patches), len(ark.ax.lines), len(ark.ax.collections)
        if flater and self.flater:
            ark.flate([self.a(f) for f in self.flater], skravur=self.skravur)
        if self.syn:
            lc = LineCollection([self.a(l) for l in self.syn], colors="k", linewidths=TYKK * PT, zorder=4,
                                capstyle="round", joinstyle="round")
            ark.ax.add_collection(lc)
        for l in self.skj:
            if len(l) > 1 and np.linalg.norm(np.diff(l, axis=0), axis=1).sum() * self.k > 0.6:
                ark.linje(self.a(l), lw=TYNN, ls=STIPLET, z=3)
        if klipp:
            sirkel = plt.Circle(klipp[0], klipp[1], transform=ark.ax.transData, fc="none", ec="none")
            ark.ax.add_patch(sirkel)
            nye = (ark.ax.patches[n0[0]:-1] + ark.ax.lines[n0[1]:] + ark.ax.collections[n0[2]:])
            for a in nye:
                a.set_clip_path(sirkel)
            ark.sirkel(klipp[0], klipp[1], lw=TYNN)

    def senterlinje(self, p0, p1):
        """Senterlinje mellom to modellpunkter."""
        self.ark.senterlinje(self.p(p0), self.p(p1))


class Sammenstilling:
    """Riss eller snitt av flere deler. deler: liste av dict(m=Manifold, skravur='45'|'-45'|None|'lager',
    tetthet=4, navn=...). Delene skjuler hverandre. I snitt skraveres hver del for seg; None gir en
    hvit snittflate (aksler, skruer og andre deler som ikke snittes på langs), 'lager' gir
    forenklet lagerfremstilling (ISO 8826-1: kryss i rektangelet)."""

    def __init__(self, ark, deler, R, U, o, k, ref=(0, 0, 0), snitt=None, vinkel=20.0, steg=None, tol=1e-6):
        self.ark, self.k = ark, k
        self.B = basis(R, U)
        self.o = np.asarray(o, float)
        self.ref = np.asarray(ref, float)
        d_snitt = None if snitt is None else float(self.B[2] @ np.asarray(snitt, float))
        steg = steg or max(0.04, 0.6 / k)
        vis = [til_vis(d["m"], self.B, d_snitt if d.get("snittes", True) else None) for d in deler]
        occ = mf.Manifold.compose([v[0] for v in vis])
        self.deler = []
        for d, v in zip(deler, vis):
            syn, skj, fl = kanter(d["m"], self.B, d.get("vinkel", vinkel), None, False, steg, tol, vis=v,
                                  okkluder=occ)
            self.deler.append(dict(d, syn=syn, flater=fl))
        self.tegn()

    def a(self, uv):
        uv = np.asarray(uv, float)
        return self.o + self.k * (uv - self.B[:2] @ self.ref)

    def p(self, x, y=None, z=None):
        q = np.asarray(x if y is None else (x, y, z), float)
        return self.o + self.k * (self.B[:2] @ (q - self.ref))

    def tegn(self):
        ark = self.ark
        for d in self.deler:
            if d["flater"]:
                sk = d.get("skravur", "45")
                polys = [self.a(f) for f in d["flater"]]
                if sk == "lager":
                    ark.flate(polys, skravur=None)
                    for pl in polys:  # forenklet lager: oppreist kryss midt i hvert snittrektangel
                        c = pl.mean(axis=0)
                        w, h = np.ptp(pl[:, 0]), np.ptp(pl[:, 1])
                        s_ = 0.3 * min(w, h)
                        ark.linje([c - [s_, 0], c + [s_, 0]], lw=TYKK, z=5)
                        ark.linje([c - [0, s_], c + [0, s_]], lw=TYKK, z=5)
                else:
                    ark.flate(polys, skravur=sk, tetthet=d.get("tetthet", 4))
            if d["syn"]:
                ark.ax.add_collection(LineCollection([self.a(l) for l in d["syn"]], colors="k",
                                                     linewidths=TYKK * PT, zorder=4, capstyle="round"))
