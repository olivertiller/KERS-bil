"""
Tegner isometriske konseptskisser av de tre drivverkskonseptene for KERS-bilen.

    Konsept 1: svinghjulet sitter direkte på drivakselen (i = 1)
    Konsept 2: to motsatt vendte kjegler med en O-ring som reim (trinnløs CVT, Evans-prinsippet)
    Konsept 3: ett fast tannhjulstrinn fra svinghjul til drivaksel (i = 4)

Bruk:
    python3 konseptskisser.py            # lagrer PDF og PNG i Images/Konsepter/

Koordinater i mm: x forover (bakkant = 0), y mot venstre side (midtlinje = 0), z opp (bakken = 0).
Alle deler er skiver eller aksler langs y, så de tegnes i rekkefølge etter y (bakerst først)
sett fra venstre side, foran og ovenfra.
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Polygon

UT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Images", "Konsepter")



def kamera(az_grader, el_grader):
    """Kamera fra venstre, foran og ovenfra. az måles fra x-aksen mot y-aksen."""
    global V, HOYRE, OPP
    az, el = math.radians(az_grader), math.radians(el_grader)
    V = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    HOYRE = np.cross(-V, [0, 0, 1])
    HOYRE /= np.linalg.norm(HOYRE)
    OPP = np.cross(HOYRE, -V)


kamera(28, 38)
LYS = np.array([0.35, 0.55, 0.85])
LYS /= np.linalg.norm(LYS)

BLEKK = "#22272b"
NAERE_HJUL_ALFA = 0.45
FARGE = {
    "stål": "#9aa6b1",
    "aksel": "#6b737a",
    "dekk": "#3a3f44",
    "felg": "#dde1e4",
    "lager": "#b9c0c6",
    "plate": "#c8bfae",
    "drivverk": "#3b78c2",
    "drivverk_lys": "#6fa2df",
    "trommel": "#d4952f",
    "gummi": "#2e3236",
    "snor": "#b5322a",
}


def proj(p):
    p = np.asarray(p, float)
    return np.array([p @ HOYRE, p @ OPP])


def skygge(farge, n):
    k = 0.62 + 0.38 * max(0.0, float(n @ LYS))
    return tuple(min(1.0, c * k) for c in to_rgb(farge))


def basis(akse):
    a = akse / np.linalg.norm(akse)
    hjelp = np.array([0, 0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(a, hjelp)
    e1 /= np.linalg.norm(e1)
    return a, e1, np.cross(a, e1)


def tannprofil(n_tenner, r_topp, r_rot):
    """Radius som funksjon av vinkel for et forenklet tannhjul."""

    def r(theta):
        t = (theta * n_tenner / (2 * math.pi)) % 1.0
        if t < 0.15:
            f = t / 0.15
        elif t < 0.45:
            f = 1.0
        elif t < 0.6:
            f = 1 - (t - 0.45) / 0.15
        else:
            f = 0.0
        return r_rot + (r_topp - r_rot) * f

    return r, max(96, n_tenner * 10)


class Del:
    """Felles for alle deler: sorteres på y-midtpunkt, deretter dybde mot kamera."""

    def nokkel(self):
        return (round(self.ymid, 3), float(self.senter @ V))


class Sylinder(Del):
    """Avkortet kjegle mellom p0 og p1. Med like radier blir det en sylinder."""

    def __init__(self, p0, p1, r0, r1=None, farge="stål", lokk=None, n=72, profil=None,
                 pynt=(), alfa=1.0):
        self.p0, self.p1 = np.asarray(p0, float), np.asarray(p1, float)
        self.r0, self.r1 = r0, r0 if r1 is None else r1
        self.farge, self.lokk = FARGE.get(farge, farge), FARGE.get(lokk, lokk) if lokk else None
        self.n, self.profil, self.pynt, self.alfa = n, profil, pynt, alfa
        if profil:
            self.n = profil[1]
        self.senter = 0.5 * (self.p0 + self.p1)
        self.ymid = self.senter[1]

    def tegn(self, ax, z):
        a, e1, e2 = basis(self.p1 - self.p0)
        th = np.linspace(0, 2 * math.pi, self.n, endpoint=False)
        if self.profil:
            rr = np.array([self.profil[0](t) for t in th])
            ring0 = self.p0 + rr[:, None] * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
            ring1 = self.p1 + rr[:, None] * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
        else:
            ret = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
            ring0, ring1 = self.p0 + self.r0 * ret, self.p1 + self.r1 * ret

        m = len(th)
        synlig = np.zeros(m, bool)
        normaler = []
        for i in range(m):
            j = (i + 1) % m
            nrm = np.cross(ring0[j] - ring0[i], ring1[i] - ring0[i])
            midt = 0.25 * (ring0[i] + ring0[j] + ring1[i] + ring1[j])
            akse_pkt = self.p0 + ((midt - self.p0) @ a) * a
            if nrm @ (midt - akse_pkt) < 0:
                nrm = -nrm
            nrm /= np.linalg.norm(nrm) + 1e-12
            normaler.append(nrm)
            synlig[i] = nrm @ V > 0

        lokk0, lokk1 = (-a) @ V > 0, a @ V > 0
        if self.alfa < 1 and not self.profil and synlig.any() and not synlig.all():
            # Gjennomsiktig: hele den synlige mantelen som én flate, ellers synes sømmene
            start = int(np.nonzero(synlig & ~np.roll(synlig, 1))[0][0])
            idx = [(start + k) % m for k in range(int(synlig.sum()) + 1)]
            poly = [proj(ring0[i]) for i in idx] + [proj(ring1[i]) for i in reversed(idx)]
            c = skygge(self.farge, normaler[idx[len(idx) // 2]])
            ax.add_patch(Polygon(poly, closed=True, fc=c, ec="none", alpha=self.alfa, zorder=z))
        else:
            for i in np.nonzero(synlig)[0]:
                j = (i + 1) % m
                poly = [proj(ring0[i]), proj(ring0[j]), proj(ring1[j]), proj(ring1[i])]
                c = skygge(self.farge, normaler[i])
                ax.add_patch(Polygon(poly, closed=True, fc=c, ec=c, lw=0.3, alpha=self.alfa, zorder=z))
        # Konturer: overgang synlig/skjult side, og kanten mot lokk som vender bort
        for i in range(m):
            j = (i + 1) % m
            if synlig[i] != synlig[i - 1]:
                self._linje(ax, [ring0[i], ring1[i]], z + 0.2)
            if synlig[i]:
                if not lokk0:
                    self._linje(ax, [ring0[i], ring0[j]], z + 0.2)
                if not lokk1:
                    self._linje(ax, [ring1[i], ring1[j]], z + 0.2)
        for vis, ring, sentrum, nrm in ((lokk0, ring0, self.p0, -a), (lokk1, ring1, self.p1, a)):
            if not vis:
                continue
            c = skygge(self.lokk or self.farge, nrm)
            ax.add_patch(Polygon([proj(p) for p in ring], closed=True, fc=c, ec=BLEKK, lw=0.7,
                                 alpha=self.alfa, zorder=z + 0.3))
            for form in self.pynt:
                pts = [sentrum + u * e1 + v * e2 for u, v in form]
                self._linje(ax, pts, z + 0.4, lw=0.45)

    def _linje(self, ax, pts, z, lw=0.7):
        xy = np.array([proj(p) for p in pts])
        ax.plot(xy[:, 0], xy[:, 1], color=BLEKK, lw=lw, alpha=max(self.alfa, 0.5), zorder=z,
                solid_capstyle="round")


class Boks(Del):
    def __init__(self, lo, hi, farge="plate", alfa=1.0):
        self.lo, self.hi = np.asarray(lo, float), np.asarray(hi, float)
        self.farge, self.alfa = FARGE.get(farge, farge), alfa
        self.senter = 0.5 * (self.lo + self.hi)
        self.ymid = self.senter[1]

    def tegn(self, ax, z):
        x0, y0, z0 = self.lo
        x1, y1, z1 = self.hi
        hj = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
              (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
        flater = [((0, 1, 2, 3), (0, 0, -1)), ((4, 5, 6, 7), (0, 0, 1)), ((0, 1, 5, 4), (0, -1, 0)),
                  ((3, 2, 6, 7), (0, 1, 0)), ((0, 3, 7, 4), (-1, 0, 0)), ((1, 2, 6, 5), (1, 0, 0))]
        for idx, nrm in flater:
            nrm = np.array(nrm, float)
            if nrm @ V <= 0:
                continue
            ax.add_patch(Polygon([proj(hj[i]) for i in idx], closed=True, fc=skygge(self.farge, nrm),
                                 ec=BLEKK, lw=0.6, alpha=self.alfa, zorder=z))


class Strek(Del):
    def __init__(self, pts, farge="snor", lw=1.2, stiplet=False, nokkel=None):
        self.pts = np.asarray(pts, float)
        self.farge, self.lw, self.stiplet = FARGE.get(farge, farge), lw, stiplet
        self.senter = self.pts.mean(axis=0)
        self.ymid = self.senter[1]
        self._nokkel = nokkel

    def nokkel(self):
        return self._nokkel or super().nokkel()

    def tegn(self, ax, z):
        xy = np.array([proj(p) for p in self.pts])
        ax.plot(xy[:, 0], xy[:, 1], color=self.farge, lw=self.lw, zorder=z,
                ls=(0, (3, 2)) if self.stiplet else "-", solid_capstyle="round")


# ---------- Byggeklosser ----------

def sirkel(r, n=48, dx=0.0, dy=0.0):
    t = np.linspace(0, 2 * math.pi, n + 1)
    return list(zip(dx + r * np.cos(t), dy + r * np.sin(t)))


def hjul(x, z, r, y_inn, y_ut, nav=1.5):
    """Printet hjul med O-ring og fem letthull, sett fra utsiden.
    Hjulene på siden nærmest kameraet tegnes halvgjennomsiktige så mekanismen synes."""
    pynt = [sirkel(r - 2.5), sirkel(nav + 2.5), sirkel(nav)]
    r_hull = 0.5 * (r - 2.5 + nav + 2.5)
    for k in range(5):
        v = 2 * math.pi * k / 5 + 0.3
        pynt.append(sirkel(0.28 * (r - nav - 5), 24, r_hull * math.cos(v), r_hull * math.sin(v)))
    return Sylinder((x, y_inn, z), (x, y_ut, z), r, farge="dekk", lokk="felg", pynt=pynt,
                    alfa=NAERE_HJUL_ALFA if y_inn > 0 else 1.0)


def aksel(x, z, r, y_fra, y_til, hull):
    """Akselbiter i mellomrommene mellom delene (hull = liste med y-intervaller som dekkes)."""
    kanter = sorted({y_fra, y_til, *[v for h in hull for v in h]})
    deler = []
    for a, b in zip(kanter[:-1], kanter[1:]):
        m = 0.5 * (a + b)
        if y_fra <= m <= y_til and not any(h[0] <= m <= h[1] for h in hull):
            deler.append(Sylinder((x, a, z), (x, b, z), r, farge="aksel", n=24))
    return deler


def lager(x, z, y0, y1):
    return Sylinder((x, y0, z), (x, y1, z), 9.5, farge="lager",
                    pynt=[sirkel(8.0), sirkel(5.0), sirkel(3.0)])


def snor(start, retning, lengde=55):
    retning = np.asarray(retning, float) / np.linalg.norm(retning)
    s = np.asarray(start, float)
    ende = s + lengde * retning
    ring = [ende + 3 * (math.cos(t) * np.array([0, 0, 1]) + math.sin(t) * retning)
            for t in np.linspace(0, 2 * math.pi, 30)]
    return [Strek([s, ende]), Strek(ring, lw=1.4)]


def sideplater(x0, x1, z0, z1, y_inn, y_ut):
    return [Boks((x0, y_inn, z0), (x1, y_ut, z1), alfa=0.32),
            Boks((x0, -y_ut, z0), (x1, -y_inn, z1), alfa=0.32)]


def reim(a, r_a, b, r_b, y, lw=2.4):
    """O-ring som reim rundt to parallelle aksler i planet y (a har størst radius).
    Bitene som ligger bak kjeglene tegnes før dem, resten etter."""
    a, b = np.array([a[0], a[1]], float), np.array([b[0], b[1]], float)
    d = b - a
    th0 = math.atan2(d[1], d[0])
    beta = math.acos((r_a - r_b) / np.linalg.norm(d))

    def pkt(c, r, v):
        return np.array([c[0] + r * math.cos(v), y, c[1] + r * math.sin(v)])

    bue_a = np.linspace(th0 + beta, th0 - beta + 2 * math.pi, 120)
    bue_b = np.linspace(th0 - beta, th0 + beta, 40)
    punkter = [(pkt(a, r_a, v), np.array([math.cos(v), 0, math.sin(v)]) @ V > 0) for v in bue_a]
    punkter += [(pkt(b, r_b, v), np.array([math.cos(v), 0, math.sin(v)]) @ V > 0) for v in bue_b]
    punkter.append(punkter[0])

    deler, biter = [], [[punkter[0]]]
    for p in punkter[1:]:
        if p[1] != biter[-1][-1][1]:
            biter.append([biter[-1][-1]])
        biter[-1].append(p)
    for bit in biter:
        foran = bit[-1][1]
        nokkel = (y + 0.01, 1e9) if foran else (y - 0.01, -1e9)
        deler.append(Strek([p for p, _ in bit], farge="gummi", lw=lw, nokkel=nokkel))
    return deler


# ---------- Konseptene ----------

def konsept1():
    """Svinghjulet sitter på drivakselen. Hjulene må være større enn svinghjulet."""
    xb, xf, zh, r = 34.0, 106.0, 34.0, 34.0
    return [
        Sylinder((xb, -3, zh), (xb, 3, zh), 30, farge="stål", pynt=[sirkel(5), sirkel(27)]),
        Sylinder((xb, -12, zh), (xb, -5, zh), 7, farge="trommel"),
        Sylinder((xb, -13.5, zh), (xb, -12, zh), 11, farge="trommel"),
        Sylinder((xb, -5, zh), (xb, -4, zh), 11, farge="trommel"),
        lager(xb, zh, 19, 25), lager(xb, zh, -25, -19),
        hjul(xb, zh, r, 26.5, 32.5, nav=3), hjul(xb, zh, r, -32.5, -26.5, nav=3),
        hjul(xf, zh, r, 26.5, 32.5), hjul(xf, zh, r, -32.5, -26.5),
        *sideplater(18, 122, 20, 50, 19, 25),
        *aksel(xb, zh, 3, -26.5, 26.5, [(-3, 3), (-13.5, -4), (19, 25), (-25, -19)]),
        *aksel(xf, zh, 1.5, -26.5, 26.5, []),
        *snor((xb, -8.5, zh + 7), (-1, 0, 0.35)),
    ]


def konsept2():
    """To motsatt vendte kjegler på parallelle aksler med en O-ring som reim.
    Når reimen flyttes langs kjeglene, endres utvekslingen trinnløst."""
    xb, xf, zh, r = 32.0, 108.0, 32.0, 32.0
    xs, zs = 65.76, 36.0
    y0, y1 = 6.0, 24.0
    rs0, rs1 = 5.0, 7.5  # kjegle på svinghjulsakselen
    rb0, rb1 = 22.0, 19.5  # kjegle på drivakselen
    yr = 15.0  # reimens posisjon i skissen

    def lerp(v0, v1):
        return v0 + (v1 - v0) * (yr - y0) / (y1 - y0)

    rs, rb = lerp(rs0, rs1), lerp(rb0, rb1)
    return [
        Sylinder((xs, -3, zs), (xs, 3, zs), 30, farge="stål", pynt=[sirkel(5), sirkel(27)]),
        Sylinder((xs, y0, zs), (xs, y1, zs), rs0, rs1, farge="drivverk", n=96),
        Sylinder((xb, y0, zh), (xb, y1, zh), rb0, rb1, farge="drivverk_lys", n=120,
                 pynt=[sirkel(4), sirkel(rb0 - 3)]),
        *reim((xb, zh), rb, (xs, zs), rs, yr),
        Sylinder((xs, -12, zs), (xs, -5, zs), 7, farge="trommel"),
        Sylinder((xs, -13.5, zs), (xs, -12, zs), 11, farge="trommel"),
        Sylinder((xs, -5, zs), (xs, -4, zs), 11, farge="trommel"),
        lager(xs, zs, 26.5, 32.5), lager(xs, zs, -32.5, -26.5),
        lager(xb, zh, 26.5, 32.5), lager(xb, zh, -32.5, -26.5),
        hjul(xb, zh, r, 33.5, 41.5), hjul(xb, zh, r, -41.5, -33.5),
        hjul(xf, zh, r, 33.5, 41.5), hjul(xf, zh, r, -41.5, -33.5),
        *sideplater(20, 118, 18, 52, 26.5, 32.5),
        *aksel(xs, zs, 3, -26.5, 26.5, [(-3, 3), (y0, y1), (-13.5, -4)]),
        *aksel(xb, zh, 1.5, -33.5, 33.5, [(y0, y1), (26.5, 32.5), (-32.5, -26.5)]),
        *aksel(xf, zh, 1.5, -33.5, 33.5, []),
        *snor((xs, -8.5, zs + 7), (-1, 0, 0.35)),
    ]


def konsept3():
    """Ett fast tannhjulstrinn, m = 0,8, 17 : 68."""
    xb, xf, zh, r = 32.0, 108.0, 32.0, 32.0
    xs, zs = 65.76, 36.0
    pinjong = tannprofil(17, 7.6, 5.8)
    tannhjul = tannprofil(68, 28.0, 26.2)
    hull_pynt = [sirkel(5.5, 24, 16 * math.cos(v), 16 * math.sin(v)) for v in np.arange(4) * math.pi / 2 + 0.4]
    return [
        Sylinder((xs, -3, zs), (xs, 3, zs), 30, farge="stål", pynt=[sirkel(5), sirkel(27)]),
        Sylinder((xs, 5, zs), (xs, 10, zs), 7.6, farge="drivverk", profil=pinjong, pynt=[sirkel(3)]),
        Sylinder((xb, 5, zh), (xb, 10, zh), 28, farge="drivverk", profil=tannhjul,
                 pynt=[sirkel(24.5), sirkel(4), *hull_pynt]),
        Sylinder((xs, -12, zs), (xs, -5, zs), 7, farge="trommel"),
        Sylinder((xs, -13.5, zs), (xs, -12, zs), 11, farge="trommel"),
        Sylinder((xs, -5, zs), (xs, -4, zs), 11, farge="trommel"),
        lager(xs, zs, 26.5, 32.5), lager(xs, zs, -32.5, -26.5),
        lager(xb, zh, 26.5, 32.5), lager(xb, zh, -32.5, -26.5),
        hjul(xb, zh, r, 33.5, 41.5), hjul(xb, zh, r, -41.5, -33.5),
        hjul(xf, zh, r, 33.5, 41.5), hjul(xf, zh, r, -41.5, -33.5),
        *sideplater(20, 118, 18, 52, 26.5, 32.5),
        *aksel(xs, zs, 3, -26.5, 26.5, [(-3, 3), (5, 10), (-13.5, -4)]),
        *aksel(xb, zh, 1.5, -33.5, 33.5, [(5, 10), (26.5, 32.5), (-32.5, -26.5)]),
        *aksel(xf, zh, 1.5, -33.5, 33.5, []),
        *snor((xs, -8.5, zs + 7), (-1, 0, 0.35)),
    ]


def omfang(deler):
    """Ytterpunktene til delene i bildeplanet."""
    pts = []
    for d in deler:
        if isinstance(d, Sylinder):
            a, e1, e2 = basis(d.p1 - d.p0)
            for p, r in ((d.p0, d.r0), (d.p1, d.r1)):
                pts += [proj(p + r * (math.cos(t) * e1 + math.sin(t) * e2)) for t in np.linspace(0, 2 * math.pi, 36)]
        elif isinstance(d, Boks):
            pts += [proj((x, y, z)) for x in (d.lo[0], d.hi[0]) for y in (d.lo[1], d.hi[1]) for z in (d.lo[2], d.hi[2])]
        else:
            pts += [proj(p) for p in d.pts]
    pts = np.array(pts)
    return pts.min(axis=0), pts.max(axis=0)


def tegn(deler, filnavn, ramme):
    # Samme ramme for alle konseptene, og figuren får rammens sideforhold,
    # så alle bildene får samme målestokk uten ekstra marg
    (x0, y0), (x1, y1) = ramme[0] - 4, ramme[1] + 4
    fig, ax = plt.subplots(figsize=(5.0 * (x1 - x0) / (y1 - y0), 5.0))
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    for z, d in enumerate(sorted(deler, key=lambda d: d.nokkel())):
        d.tegn(ax, z)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    os.makedirs(UT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(UT, f"{filnavn}.{ext}"), dpi=300, facecolor="white")
    plt.close(fig)


def main():
    konsepter = {"konsept1_direkte": konsept1(), "konsept2_cvt": konsept2(),
                 "konsept3_tannhjul": konsept3()}
    grenser = [omfang(d) for d in konsepter.values()]
    ramme = (np.min([g[0] for g in grenser], axis=0), np.max([g[1] for g in grenser], axis=0))
    for navn, deler in konsepter.items():
        tegn(deler, navn, ramme)
    print(f"Lagret i {os.path.normpath(UT)}")


if __name__ == "__main__":
    main()
