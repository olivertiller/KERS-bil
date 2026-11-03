"""Laserkuttet MDF-kropp v1 for KERS-bilen (Epilog Fusion M2 via FlexiDESIGNER).

Lager kuttfiler (SVG, DXF og PDF i målestokk 1:1, mm) og en forhåndsvisning:
    python3 laserkutt_kropp.py                    # standardverdier
    python3 laserkutt_kropp.py --t 3.1 --lager 18.9 --spalte 2.95

Koordinater som i drivverk-tegningen: x fra bakkant, y fra midtlinjen
(+y = tannhjulssiden, -y = trommelsiden), z fra bakken. Alle mål i mm.

Rød strek = kutt, blå strek = gravering (tekst og markeringer, vektor med lav effekt).

Kroppen:
  * 2 innervegger + 2 yttervegger, limes sammen to og to til 6 mm vegg (= bredden på 626-lageret).
    Hver vegg har tre lagerhull på z = 34: x = 34 (bakaksel), 70 (svinghjul), 106 (foraksel).
  * Bunnplate med åpning for svinghjul og tannhjul. Bunnplata ligger mellom veggene.
  * Frontskott og bakskott med tapper inn i innerveggene og bunnplata.
  * Holder for Hall-sensoren ved trommelens ytre flens.

Alle tapper går sideveis inn i veggene, så en vegg kan skyves av og på langs akslene.
Hver vegg holdes av to M3 x 12 inn i T-spor i kanten av bunnplata. Skruene sitter nederst
mellom bak- og forhjulet, det eneste stedet skruehodene går klar av hjulene. Den ene veggen
kan limes fast, den andre skrus på etter at svinghjulet er satt inn.
"""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle as MplCircle
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath

UT = Path(__file__).resolve().parent.parent / "Kers-bil" / "Laserkutt"

# ---------------------------------------------------------------- parametere
P = dict(
    t=3.0,  # platetykkelse (mål med skyvelære!)
    lager=18.8,  # tegnet lagerhull; laseren brenner bort ca. 0,15 → ca. 18,95 ferdig
    spalte=2.85,  # tegnet spaltebredde for tapper (ca. t - kerf)
    x_akser=(34.0, 70.0, 106.0),
    z_aksel=34.0,
    # vegger (sideriss, x og z)
    x0=2.0,
    x1=132.0,
    z_vegg=5.0,  # underkant vegger
    z_bunn=8.0,  # underkant bunnplate
    z_topp=50.0,  # overkant vegger og frontskott
    y_inn=26.5,  # innerflate vegg (= halv innvendig bredde = halv bredde bunnplate)
    # skott
    x_front=125.0,  # bakre flate frontskott
    x_bak=6.0,  # bakre flate bakskott
    z_bak_topp=30.0,
    # åpning i bunnplata
    aapning=(10.0, 95.0, -6.0, 13.0),  # x0, x1, y0, y1
    # veggskruer: M3 x 12 gjennom veggen (6) og inn i T-spor i kanten av bunnplata
    x_skruer=(64.0, 76.0),
    m3_hull=3.2,
    mutter_b=5.5,  # tegnet lomme for M3-mutter (5,5 over flatene)
    mutter_h=2.5,  # tegnet lommebredde (mutter 2,4 tykk)
    mutter_inn=4.0,  # mutterens senter fra platekanten
    tspor_dybde=7.5,  # skruesporets dybde fra kanten (skruen rekker 6 inn)
    skruehode=5.6,  # Ø skruehode (linsehode M3)
    # sensorholder
    y_sensor=-17.0,  # flaten sensoren limes på (trommelens ytre flens er på y = -13,5)
)


# ---------------------------------------------------------------- geometri
class Del:
    """En del: lukkede polylinjer og sirkler i delens egne koordinater (u, v)."""

    def __init__(self, navn, antall=1):
        self.navn, self.antall = navn, antall
        self.kutt = []  # lukkede polylinjer [(u, v), ...]
        self.sirkler = []  # (u, v, d)
        self.grav = []  # polylinjer (tekst og merker), åpne eller lukkede

    def bbox(self):
        pts = [p for pl in self.kutt for p in pl]
        pts += [(u - d / 2, v - d / 2) for u, v, d in self.sirkler]
        pts += [(u + d / 2, v + d / 2) for u, v, d in self.sirkler]
        us, vs = zip(*pts)
        return min(us), min(vs), max(us), max(vs)


def rekt(u0, v0, u1, v1):
    return [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]


def spalte_u(uc, v0, v1, b):
    """Spalte med bredde b langs u, sentrert på uc, fra v0 til v1."""
    return rekt(uc - b / 2, v0, uc + b / 2, v1)


def spalte_v(vc, u0, u1, b):
    return rekt(u0, vc - b / 2, u1, vc + b / 2)


def bue(cu, cv, r, a0, a1, n=8):
    return [(cu + r * math.cos(math.radians(a0 + (a1 - a0) * i / n)),
             cv + r * math.sin(math.radians(a0 + (a1 - a0) * i / n))) for i in range(n + 1)]


def kant_med_tapper(p0, p1, tapper, dybde, normal):
    """Punkter fra p0 mot p1 (p1 ikke med). tapper = [(s0, s1)] avstand langs kanten."""
    (x0, y0), (x1, y1) = p0, p1
    L = math.hypot(x1 - x0, y1 - y0)
    ex, ey = (x1 - x0) / L, (y1 - y0) / L
    nx, ny = normal
    pts = [p0]
    for s0, s1 in sorted(tapper):
        a = (x0 + ex * s0, y0 + ey * s0)
        b = (x0 + ex * s1, y0 + ey * s1)
        pts += [a, (a[0] + nx * dybde, a[1] + ny * dybde), (b[0] + nx * dybde, b[1] + ny * dybde), b]
    return pts


def kontur(u0, v0, u1, v1, r_topp=0.0, bunn=(), topp=(), venstre=(), hoyre=(), d=3.0):
    """Rektangel med tapper på kantene og avrundede øvre hjørner (radius r_topp).
    Tapper oppgis i delens koordinater: bunn/topp som (u0, u1), venstre/høyre som (v0, v1)."""
    pts = kant_med_tapper((u0, v0), (u1, v0), [(a - u0, b - u0) for a, b in bunn], d, (0, -1))
    pts += kant_med_tapper((u1, v0), (u1, v1 - r_topp), [(a - v0, b - v0) for a, b in hoyre], d, (1, 0))
    if r_topp:
        pts += bue(u1 - r_topp, v1 - r_topp, r_topp, 0, 90)[:-1]
    pts += kant_med_tapper((u1 - r_topp, v1), (u0 + r_topp, v1), [(u1 - r_topp - b, u1 - r_topp - a) for a, b in topp], d, (0, 1))
    if r_topp:
        pts += bue(u0 + r_topp, v1 - r_topp, r_topp, 90, 180)[:-1]
    pts += kant_med_tapper((u0, v1 - r_topp), (u0, v0), [(v1 - r_topp - b, v1 - r_topp - a) for a, b in venstre], d, (-1, 0))
    return pts


def kant_med_tspor(p0, p1, tapper, tspor, dybde, normal, p):
    """Som kant_med_tapper, men med T-spor (skruespor + mutterlomme) innover fra kanten
    i posisjonene tspor (avstand langs kanten)."""
    (x0, y0), (x1, y1) = p0, p1
    L = math.hypot(x1 - x0, y1 - y0)
    ex, ey = (x1 - x0) / L, (y1 - y0) / L
    nx, ny = normal
    c, n = p["m3_hull"] / 2, p["mutter_b"] / 2
    m0, m1 = p["mutter_inn"] - p["mutter_h"] / 2, p["mutter_inn"] + p["mutter_h"] / 2
    dc = p["tspor_dybde"]

    def P(s, h):  # h > 0 er utover
        return (x0 + ex * s + nx * h, y0 + ey * s + ny * h)

    ting = [(s0, "tapp", s1) for s0, s1 in tapper] + [(s - n, "tspor", s) for s in tspor]
    pts = [p0]
    for _, typ, s in sorted(ting):
        if typ == "tapp":
            s0 = _
            pts += [P(s0, 0), P(s0, dybde), P(s, dybde), P(s, 0)]
        else:
            pts += [P(s - c, 0), P(s - c, -m0), P(s - n, -m0), P(s - n, -m1), P(s - c, -m1), P(s - c, -dc),
                    P(s + c, -dc), P(s + c, -m1), P(s + n, -m1), P(s + n, -m0), P(s + c, -m0), P(s + c, 0)]
    return pts


def avrundet_rekt(u0, v0, u1, v1, r):
    pts = bue(u1 - r, v0 + r, r, -90, 0)[:-1] + bue(u1 - r, v1 - r, r, 0, 90)[:-1]
    pts += bue(u0 + r, v1 - r, r, 90, 180)[:-1] + bue(u0 + r, v0 + r, r, 180, 270)[:-1]
    return pts


_FONT = FontProperties(family="DejaVu Sans", weight="bold")


def tekst(s, u, v, h=3.5, midt=True):
    """Tekst som omriss (polylinjer), slik at den overlever import uten fonter."""
    tp = TextPath((0, 0), s, size=h, prop=_FONT)
    bb = tp.get_extents()
    du = u - (bb.x0 + bb.width / 2 if midt else bb.x0)
    dv = v - (bb.y0 + bb.height / 2)
    return [[(x + du, y + dv) for x, y in poly] for poly in tp.to_polygons(closed_only=True)]


def pil(u, v, lengde=14, retning=1):
    """Pil langs +u (retning=1) eller -u."""
    s = retning
    return [[(u - s * lengde / 2, v), (u + s * lengde / 2, v)],
            [(u + s * (lengde / 2 - 3), v + 1.8), (u + s * lengde / 2, v), (u + s * (lengde / 2 - 3), v - 1.8)]]


# ---------------------------------------------------------------- delene
def lag_deler(p):
    t, sp, dl = p["t"], p["spalte"], p["lager"]
    z_plate = p["z_bunn"] + t  # overkant bunnplate = underkant vegger og skott
    x0, x1, zt, zb = p["x0"], p["x1"], p["z_topp"], p["z_bunn"]
    yi = p["y_inn"]

    # tapper fra bunnplata inn i innerveggene (x-intervaller); ikke ved skruene på x 64 og 76
    vegg_tapper = [(14.0, 30.0), (38.0, 54.0), (86.0, 102.0), (108.0, 120.0)]
    z_skrue = p["z_bunn"] + t / 2  # midt i bunnplata
    skruehull = [(x, z_skrue, p["m3_hull"]) for x in p["x_skruer"]]
    # bakskott og frontskott: tapper inn i veggene (z-intervaller) og bunnplata (y-intervaller)
    xf = p["x_front"] + t / 2
    xb = p["x_bak"] + t / 2
    front_side = (22.0, 38.0)
    bak_side = (14.0, 24.0)
    front_bunn = [(-20.0, -8.0), (8.0, 20.0)]
    bak_bunn = [(-22.0, -10.0), (16.0, 24.0)]
    sensor_x = (60.0, 80.0)
    sensor_z = 31.0
    sensor_bunn = [(68.0, 72.0)]  # mellom T-sporene for veggskruene
    ys = p["y_sensor"] - t / 2  # midt i sensorholderen

    deler = []

    # Innervegg (sideriss: u = x, v = z)
    d = Del("Innervegg", 2)
    d.kutt.append(kontur(x0, p["z_vegg"], x1, zt, r_topp=6))
    d.kutt.append(spalte_u(xf, *front_side, sp))
    d.kutt.append(spalte_u(xb, *bak_side, sp))
    for a, b in vegg_tapper:
        d.kutt.append(spalte_v(z_skrue, a, b, sp))
    d.sirkler += [(x, p["z_aksel"], dl) for x in p["x_akser"]] + skruehull
    d.grav += tekst("INNERVEGG", 52, 46, 3.0) + tekst("2 stk", 88, 46, 3.0)
    d.grav += pil(116, 18) + tekst("FORAN", 116, 23, 2.8)
    deler.append(d)

    # Yttervegg
    d = Del("Yttervegg", 2)
    d.kutt.append(kontur(x0, p["z_vegg"], x1, zt, r_topp=6))
    d.sirkler += [(x, p["z_aksel"], dl) for x in p["x_akser"]] + skruehull
    d.grav += tekst("YTTERVEGG", 52, 46, 3.0) + tekst("2 stk", 88, 46, 3.0)
    deler.append(d)

    # Bunnplate (planriss: u = x, v = y)
    d = Del("Bunnplate")
    xs_ = p["x_skruer"]
    pts = kant_med_tspor((x0, -yi), (x1, -yi), [(a - x0, b - x0) for a, b in vegg_tapper],
                         [x - x0 for x in xs_], t, (0, -1), p)
    pts += [(x1, -yi)]
    pts += kant_med_tspor((x1, yi), (x0, yi), [(x1 - b, x1 - a) for a, b in vegg_tapper],
                          [x1 - x for x in xs_], t, (0, 1), p)
    pts += [(x0, yi)]
    d.kutt.append(pts)
    ax0, ax1, ay0, ay1 = p["aapning"]
    d.kutt.append(avrundet_rekt(ax0, ay0, ax1, ay1, 2))
    for a, b in front_bunn:
        d.kutt.append(spalte_u(xf, a, b, sp))
    for a, b in bak_bunn:
        d.kutt.append(spalte_u(xb, a, b, sp))
    for a, b in sensor_bunn:
        d.kutt.append(spalte_v(ys, a, b, sp))
    d.grav += tekst("BUNN (oppside)", 106, 22, 3.0)
    d.grav += pil(113, 6) + tekst("FORAN", 113, 0, 2.8)
    d.grav += tekst("tannhjul", 52, 18.5, 2.6) + tekst("trommel", 30, -12, 2.6) + tekst("sensor", 70, -14, 2.6)
    deler.append(d)

    # Frontskott (sett bakfra: u = y, v = z)
    d = Del("Frontskott")
    d.kutt.append(kontur(-yi, z_plate, yi, zt, r_topp=4, bunn=front_bunn,
                         venstre=[front_side], hoyre=[front_side], d=t))
    d.grav += tekst("FRONT", 0, 40, 3.2)
    deler.append(d)

    # Bakskott (sett bakfra: u = y, v = z) med hakk for tannhjulet på y = +2 til +13
    d = Del("Bakskott")
    zbt = p["z_bak_topp"]
    hakk = (2.0, 13.0, 15.0)  # y0, y1, z-bunn
    pts = kant_med_tapper((-yi, z_plate), (yi, z_plate), [(a + yi, b + yi) for a, b in bak_bunn], t, (0, -1))
    pts += kant_med_tapper((yi, z_plate), (yi, zbt), [(bak_side[0] - z_plate, bak_side[1] - z_plate)], t, (1, 0))
    pts += [(yi, zbt), (hakk[1], zbt), (hakk[1], hakk[2]), (hakk[0], hakk[2]), (hakk[0], zbt)]
    pts += kant_med_tapper((-yi, zbt), (-yi, z_plate), [(zbt - bak_side[1], zbt - bak_side[0])], t, (-1, 0))
    d.kutt.append(pts)
    d.grav += tekst("BAK", -12, 22, 3.2) + tekst("hakk = tannhjul", -6, 16.5, 2.4)
    deler.append(d)

    # Sensorholder (sideriss: u = x, v = z). Sensoren limes med merket front mot trommelen.
    d = Del("Sensorholder")
    d.kutt.append(kontur(sensor_x[0], z_plate, sensor_x[1], sensor_z, r_topp=2, bunn=sensor_bunn, d=t))
    zm = p["z_aksel"] - 7  # magnetene sitter på radius 7
    d.grav.append(rekt(70 - 2.05, zm - 1.5, 70 + 2.05, zm + 1.5) + [(70 - 2.05, zm - 1.5)])
    d.grav += tekst("HALL", 70, 17, 2.6)
    deler.append(d)

    return deler


def lag_lagerplate():
    """Liten plate med lagerhull Ø18,5–19,2 (tegnet mål), to rader à fire.
    Graveringen viser tegnet diameter og markerer oppsiden (laseren lager hullet litt konisk)."""
    d = Del("Lagerplate")
    verdier = [18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2]
    d.kutt.append(avrundet_rekt(0, 0, 97, 59, 3))
    for i, v in enumerate(verdier):
        rad, kol = divmod(i, 4)
        u = 12 + 24 * kol
        vh, vt = (46.0, 32.5) if rad == 0 else (17.0, 3.5)
        d.sirkler.append((u, vh, v))
        d.grav += tekst(f"{v:.1f}", u, vt, 3.0)
    return [d]


def lag_testdeler(p):
    t = p["t"]
    deler = []

    d = Del("Lagertest")
    # Se også lag_lagerplate() for en egen plate med flere diametre.
    verdier = [18.6, 18.7, 18.8, 18.9, 19.0]
    d.kutt.append(avrundet_rekt(0, 0, 22 * len(verdier) + 4, 31, 2))
    for i, v in enumerate(verdier):
        u = 13 + 22 * i
        d.sirkler.append((u, 18, v))
        d.grav += tekst(f"{v:.1f}", u, 4.5, 3.0)
    deler.append(d)

    d = Del("Spaltetest")
    verdier = [2.7, 2.8, 2.9, 3.0, 3.1]
    d.kutt.append(avrundet_rekt(0, 0, 14 * len(verdier) + 4, 24, 2))
    for i, v in enumerate(verdier):
        u = 9 + 14 * i
        d.kutt.append(spalte_u(u, 9, 21, v))
        d.grav += tekst(f"{v:.1f}", u, 4.5, 2.8)
    deler.append(d)

    d = Del("Tunge")  # stikkes i spaltene; har samme tykkelse som platene
    d.kutt.append(rekt(0, 0, 40, 10))
    deler.append(d)

    for a in (36.0, 36.1, 36.2):
        d = Del(f"Akseavstand {a:.1f}")
        d.kutt.append(avrundet_rekt(0, 0, a + 26, 30, 3))
        d.sirkler += [(13, 17, p["lager"]), (13 + a, 17, p["lager"])]
        d.grav += tekst(f"a = {a:.1f}", 13 + a / 2, 4, 2.8)
        deler.append(d)
    return deler


# ---------------------------------------------------------------- utlegg og eksport
def legg_ut(deler, bredde, luft=4.0):
    """Enkel hylle-pakking. Returnerer [(del, du, dv)] og arkstørrelse."""
    plasser, u, v, rad_h, maks_u = [], luft, luft, 0.0, 0.0
    for d in deler:
        for _ in range(d.antall):
            b0, b1, b2, b3 = d.bbox()
            w, h = b2 - b0, b3 - b1
            if u + w + luft > bredde and u > luft:
                u, v, rad_h = luft, v + rad_h + luft, 0.0
            plasser.append((d, u - b0, v - b1))
            u += w + luft
            rad_h = max(rad_h, h)
            maks_u = max(maks_u, u)
    return plasser, maks_u, v + rad_h + luft


def eksporter(plasser, W, H, stamme):
    # SVG-koordinater har y nedover, så v speiles.
    def fy(v):
        return H - v

    kutt, grav = [], []
    for d, du, dv in plasser:
        for pl in d.kutt:
            kutt.append(("pl", [(x + du, y + dv) for x, y in pl], True))
        for x, y, dia in d.sirkler:
            kutt.append(("c", (x + du, y + dv, dia), True))
        for pl in d.grav:  # graveringer er allerede lukket der de skal være det
            grav.append(("pl", [(x + du, y + dv) for x, y in pl], False))

    # SVG
    def svg_el(el, farge):
        typ, data, lukket = el
        stil = f'fill="none" stroke="{farge}" stroke-width="0.01"'
        if typ == "c":
            x, y, dia = data
            return f'<circle cx="{x:.3f}" cy="{fy(y):.3f}" r="{dia / 2:.3f}" {stil}/>'
        pts = " ".join(f"{x:.3f},{fy(y):.3f}" for x, y in data)
        return f'<{"polygon" if lukket else "polyline"} points="{pts}" {stil}/>'

    svg = [f'<?xml version="1.0" encoding="UTF-8"?>',
           f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.2f}mm" height="{H:.2f}mm" viewBox="0 0 {W:.2f} {H:.2f}">',
           '<g id="KUTT">'] + [svg_el(e, "#FF0000") for e in kutt] + ['</g>', '<g id="GRAVER">']
    svg += [svg_el(e, "#0000FF") for e in grav] + ['</g>', '</svg>']
    (UT / f"{stamme}.svg").write_text("\n".join(svg), encoding="utf-8")

    # DXF R12 (mm). Lag KUTT (rød, farge 1) og GRAVER (blå, farge 5).
    dx = ["0", "SECTION", "2", "HEADER", "9", "$INSUNITS", "70", "4", "9", "$MEASUREMENT", "70", "1",
          "0", "ENDSEC", "0", "SECTION", "2", "TABLES", "0", "TABLE", "2", "LAYER", "70", "2"]
    for navn, farge in (("KUTT", 1), ("GRAVER", 5)):
        dx += ["0", "LAYER", "2", navn, "70", "0", "62", str(farge), "6", "CONTINUOUS"]
    dx += ["0", "ENDTAB", "0", "ENDSEC", "0", "SECTION", "2", "ENTITIES"]
    for lag, farge, elementer in (("KUTT", 1, kutt), ("GRAVER", 5, grav)):
        for typ, data, lukket in elementer:
            if typ == "c":
                x, y, dia = data
                dx += ["0", "CIRCLE", "8", lag, "62", str(farge), "10", f"{x:.4f}", "20", f"{y:.4f}", "30", "0",
                       "40", f"{dia / 2:.4f}"]
            else:
                dx += ["0", "POLYLINE", "8", lag, "62", str(farge), "66", "1", "70", "1" if lukket else "0"]
                for x, y in data:
                    dx += ["0", "VERTEX", "8", lag, "10", f"{x:.4f}", "20", f"{y:.4f}", "30", "0"]
                dx += ["0", "SEQEND", "8", lag]
    dx += ["0", "ENDSEC", "0", "EOF"]
    (UT / f"{stamme}.dxf").write_text("\n".join(dx) + "\n", encoding="ascii")

    # PDF i 1:1 med hårstreker
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    for elementer, farge in ((kutt, "#FF0000"), (grav, "#0000FF")):
        for typ, data, lukket in elementer:
            if typ == "c":
                x, y, dia = data
                ax.add_patch(MplCircle((x, y), dia / 2, fill=False, ec=farge, lw=0.05))
            else:
                xs, ys = zip(*(list(data) + ([data[0]] if lukket else [])))
                ax.plot(xs, ys, color=farge, lw=0.05)
    fig.savefig(UT / f"{stamme}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- kontroll og forhåndsvisning
def kontroll(p):
    """Regner klaringer mot drivverket fra tegningsgrunnlaget."""
    t, z_a = p["t"], p["z_aksel"]
    zpl = p["z_bunn"] + t
    r_sv, r_tann, r_hjul = 30.0, 30.0, 34.0
    x_bak, x_sv, x_for = p["x_akser"]
    ax0, ax1, ay0, ay1 = p["aapning"]
    ut = {}
    # bunnplata (z_bunn..zpl) mot svinghjul (y ±3) og tannhjul (y +5..+10)
    h = z_a - zpl
    halv = math.sqrt(r_sv**2 - h**2)
    ut["åpning foran svinghjul (x)"] = ax1 - (x_sv + halv)
    ut["åpning bak tannhjul (x)"] = (x_bak - math.sqrt(r_tann**2 - h**2)) - ax0
    ut["åpning sideveis svinghjul (y)"] = -3.0 - ay0
    ut["åpning sideveis tannhjul (y)"] = ay1 - 10.0
    # bakskott: bakre del av tannhjulet passerer hakket (z >= 15) ved x = x_bak + t
    xs = p["x_bak"] + t
    ut["tannhjul over hakkbunn i bakskott (z)"] = (z_a - math.sqrt(r_tann**2 - (x_bak - xs) ** 2)) - 15.0
    ut["hakk sideveis mot tannhjul (y)"] = min(5.0 - 2.0, 13.0 - 10.0)
    ut["frontskott bak svinghjul (x)"] = p["x_front"] - (x_sv + r_sv)
    ut["sensorholder mot ytre trommelflens (y)"] = -13.5 - p["y_sensor"]
    ut["luft sensorflate–flens med 1,5 mm sensor"] = -13.5 - (p["y_sensor"] + 1.5)
    zs = p["z_bunn"] + t / 2
    ut["skruehode mot bakhjul (sideriss)"] = math.hypot(p["x_skruer"][0] - x_bak, zs - z_a) - r_hjul - p["skruehode"] / 2
    ut["skruehode mot forhjul (sideriss)"] = math.hypot(x_for - p["x_skruer"][1], zs - z_a) - r_hjul - p["skruehode"] / 2
    ut["hjul mot yttervegg"] = 33.5 - (p["y_inn"] + 2 * t)
    ut["vegg ytterflate mot svinghjulsakselens ende (32,5)"] = 32.5 - (p["y_inn"] + 2 * t)
    ut["bakkeklaring bunnplate / vegger"] = min(p["z_bunn"], p["z_vegg"])
    ut["høyde vegg (maks 70)"] = p["z_topp"]
    ut["bredde over hjul (maks 85)"] = 2 * (p["y_inn"] + 2 * t + 1 + 8)
    return ut


def forhandsvisning(p, plasser, W, H):
    t = p["t"]
    zpl = p["z_bunn"] + t
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(9, 15), gridspec_kw=dict(height_ratios=[1, 1.05, H / W * 1.6]))
    blaa, oransje, gra = "#1f5fa0", "#b4541a", "#888"

    # Sideriss
    a1.set_title("Sideriss (fra høyre), kropp i brunt, drivverk i blått")
    a1.add_patch(Rectangle((0, 0), 140, 70, fill=False, ls="--", ec=oransje, lw=0.8))
    a1.axhline(0, color="k", lw=1)
    a1.add_patch(Rectangle((p["x0"], p["z_vegg"]), p["x1"] - p["x0"], p["z_topp"] - p["z_vegg"], fc="#d9b38c", ec="#7a5230", alpha=.6))
    for x in p["x_skruer"]:
        a1.add_patch(MplCircle((x, p["z_bunn"] + t / 2), p["skruehode"] / 2, fc="#333"))
    a1.add_patch(Rectangle((p["x0"], p["z_bunn"]), p["x1"] - p["x0"], t, fc="#c49a6c", ec="#7a5230"))
    for x in p["x_akser"]:
        a1.add_patch(MplCircle((x, 34), 34, fill=False, ec="k", lw=1))
        a1.add_patch(MplCircle((x, 34), 9.5, fill=False, ec="k", ls=":"))
        a1.add_patch(MplCircle((x, 34), 3, fc="k"))
    a1.add_patch(MplCircle((70, 34), 30, fc=gra, alpha=.5, ec="k"))
    a1.add_patch(MplCircle((34, 34), 30, fill=False, ec=blaa, lw=1.2))
    a1.add_patch(MplCircle((70, 34), 7.6, fc=blaa, alpha=.5))
    a1.add_patch(Rectangle((p["x_front"], zpl), t, p["z_topp"] - zpl, fc="#7a5230"))
    a1.add_patch(Rectangle((p["x_bak"], zpl), t, p["z_bak_topp"] - zpl, fc="#7a5230"))
    a1.text(70, 6.2, "veggskruer", fontsize=7, ha="center", va="top")
    a1.text(p["x_bak"] + 4, 26, "bakskott\n(hakk for\ntannhjul)", fontsize=7)
    a1.set_xlim(-8, 148); a1.set_ylim(-4, 74); a1.set_aspect("equal")

    # Planriss
    a2.set_title("Planriss (ovenfra, foran til høyre)")
    a2.add_patch(Rectangle((0, -42.5), 140, 85, fill=False, ls="--", ec=oransje, lw=0.8))
    yb = p["y_inn"]
    a2.add_patch(Rectangle((p["x0"], -yb), p["x1"] - p["x0"], 2 * yb, fc="#e8d3bb", ec="#7a5230"))
    ax0, ax1, ay0, ay1 = p["aapning"]
    a2.add_patch(Rectangle((ax0, ay0), ax1 - ax0, ay1 - ay0, fc="white", ec="#7a5230"))
    for s in (1, -1):
        y0 = p["y_inn"] if s > 0 else -p["y_inn"] - 2 * t
        a2.add_patch(Rectangle((p["x0"], y0), p["x1"] - p["x0"], 2 * t, fc="#c49a6c", ec="#7a5230"))
        for x in p["x_akser"]:
            a2.add_patch(Rectangle((x - 9.5, y0), 19, 6, fc="white", ec="k", lw=.6))
        for x in (p["x_akser"][0], p["x_akser"][2]):
            a2.add_patch(Rectangle((x - 34, 33.5 if s > 0 else -41.5), 68, 8, fc="#555", ec="k"))
        for x in p["x_skruer"]:
            a2.add_patch(Rectangle((x - 1.5, (p["y_inn"] + 2 * t) if s > 0 else -p["y_inn"] - 2 * t - 2), 3, 2, fc="#333"))
    for x, (y0, y1) in zip(p["x_akser"], ((-41.5, 41.5), (-32.5, 32.5), (-41.5, 41.5))):
        a2.add_patch(Rectangle((x - 3, y0), 6, y1 - y0, fc="#bbb", ec="k", lw=.5))
    a2.add_patch(Rectangle((40, -3), 60, 6, fc=gra, ec="k"))
    a2.add_patch(Rectangle((4, 5), 60, 5, fc=blaa, alpha=.6))
    a2.add_patch(Rectangle((62.4, 5), 15.2, 5, fc=blaa, alpha=.6))
    a2.add_patch(Rectangle((27, 10), 14, 6, fc=blaa, alpha=.35))
    a2.add_patch(Rectangle((63, 10), 14, 6, fc=blaa, alpha=.35))
    a2.add_patch(Rectangle((59, -13.5), 22, 9.5, fc="#ddd", ec="k"))
    a2.add_patch(Rectangle((60, p["y_sensor"] - t), 20, t, fc="#7a5230"))
    a2.add_patch(Rectangle((p["x_front"], -p["y_inn"]), t, 2 * p["y_inn"], fc="#7a5230"))
    a2.add_patch(Rectangle((p["x_bak"], -p["y_inn"]), t, 2 * p["y_inn"], fc="#7a5230"))
    a2.text(72, -22, "sensorholder", fontsize=7, ha="center")
    a2.text(30, 20, "+y tannhjulsside", fontsize=7); a2.text(14, -26, "-y trommelside", fontsize=7)
    a2.set_xlim(-8, 148); a2.set_ylim(-46, 46); a2.set_aspect("equal")

    # Kuttark
    a3.set_title(f"Kuttark {W:.0f} × {H:.0f} mm (rødt = kutt, blått = gravering)")
    for d, du, dv in plasser:
        for pl in d.kutt:
            xs, ys = zip(*(pl + [pl[0]]))
            a3.plot([x + du for x in xs], [y + dv for y in ys], "r", lw=.6)
        for x, y, dia in d.sirkler:
            a3.add_patch(MplCircle((x + du, y + dv), dia / 2, fill=False, ec="r", lw=.6))
        for pl in d.grav:
            xs, ys = zip(*pl)
            a3.plot([x + du for x in xs], [y + dv for y in ys], "b", lw=.3)
    a3.add_patch(Rectangle((0, 0), W, H, fill=False, ec="k", lw=.5))
    a3.set_xlim(-3, W + 3); a3.set_ylim(-3, H + 3); a3.set_aspect("equal")
    for a in (a1, a2, a3):
        a.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(UT / "forhandsvisning_kropp_v1.png", dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t", type=float, default=P["t"], help="målt platetykkelse [mm]")
    ap.add_argument("--lager", type=float, default=P["lager"], help="tegnet lagerhull fra lagertesten [mm]")
    ap.add_argument("--spalte", type=float, default=P["spalte"], help="tegnet spaltebredde fra spaltetesten [mm]")
    a = ap.parse_args()
    p = dict(P, t=a.t, lager=a.lager, spalte=a.spalte)
    UT.mkdir(parents=True, exist_ok=True)

    deler = lag_deler(p)
    plasser, W, H = legg_ut(deler, bredde=290)
    eksporter(plasser, W, H, "kropp_v1")
    forhandsvisning(p, plasser, W, H)
    print(f"kropp_v1: {W:.0f} × {H:.0f} mm, {sum(d.antall for d in deler)} deler")

    tdeler = lag_testdeler(p)
    tpl, tW, tH = legg_ut(tdeler, bredde=200)
    eksporter(tpl, tW, tH, "test_v1")
    print(f"test_v1:  {tW:.0f} × {tH:.0f} mm, {len(tdeler)} deler")

    lpl, lW, lH = legg_ut(lag_lagerplate(), bredde=200)
    eksporter(lpl, lW, lH, "lagerplate")
    print(f"lagerplate: {lW:.0f} × {lH:.0f} mm")

    print("\nKlaringer [mm]:")
    for k, v in kontroll(p).items():
        print(f"  {k:<52} {v:6.1f}")
    print(f"\nFiler i {UT}")


if __name__ == "__main__":
    main()
