"""Snortromler TR-1 til TR-3 for KERS-bilen, som 3MF-filer (mm) klare for slicing.

    python3 Code/snortrommel.py                   # TR-1 Ø12, TR-2 Ø14, TR-3 Ø18
    python3 Code/snortrommel.py --kjerner 13 15   # andre kjernediametre
    python3 Code/snortrommel.py --boring 6.3      # hvis boringen blir for trang

Krever manifold3d:  python3 -m pip install manifold3d  (bruk gjerne en venv).

Trommelen printes stående med indre flens mot platen. Koordinater for delen: z langs aksen,
z = 0 er flaten som ligger mot svinghjulet. I bilen tilsvarer det y = -3 - z (trommelsiden).

Fra platen og opp:
  * Indre flens, 3,5 mm tykk, med utsparing for sveisekransen rundt akselen på svinghjulet
    (Ø13, maks 3 mm høy). Utsparingen er Ø14,5 rett inn 1 mm, så 45° innover mot boringen
    (5,1 mm dyp ved boringen). 45° gjør at den printes uten support med flensen mot platen.
  * Kjerne, 6 mm bred. Her ligger snora i ett lag (20 cm = 3,4-4,9 vindinger).
  * Ytre flens med 45° skrå innside, så den kan printes uten support.
  * Klemnav Ø20 x 7,5 med to slisser og to M3-skruer med mutter på tvers av slissene.
    Skruene sitter 180° fra hverandre, slik at trommelen er balansert. Ingen frest flate trengs.

Snora festes med en løkke på en krok som stikker ut fra kjernen (som en knaggkrok på en vegg):
  * En radiell stilk med en finger på toppen som peker MOT rotasjonsretningen under opptrekk.
  * Løkka legges over kroken og ligger rundt stilkfoten, under fingeren. Snora går videre
    i rotasjonsretningen og vikles over kroken.
  * Under opplading og opptrekk trekker snora løkka i rotasjonsretningen, så den trykkes mot
    stilken og fingeren hindrer den i å gli ut. Når snora er helt avviklet og trommelen spinner
    videre, snur trekket seg bakover og utover, og løkka glir av over fingertuppen av seg selv.
  * Kroken står midt i snorbanen med klaring til begge flensene, slik at løkka kommer rundt
    stilken. Undersiden skrår 45°, så den printes uten support.
Rotasjonen under opptrekk er mot klokka sett fra navenden (snora går ut fra toppen av
trommelen og bakover i bilen).
"""

import argparse
import math
import zipfile
from pathlib import Path

import manifold3d as mf
import numpy as np

UT = Path(__file__).resolve().parent.parent / "Bygg" / "Snortrommel"
SEG = 160  # segmenter i sirkler
RHO_PLA = 1.24e-3  # g/mm^3

P = dict(
    boring=6.2,  # Ø boring, samme som de andre navene
    d_utsparing=14.5,  # Ø utsparing for sveisekransen (Ø13 x maks 3 mm) ved svinghjulet
    rett_utsparing=1.0,  # rett del av utsparingen før 45°-taket
    t_indre=3.5,  # indre flens (tykk nok til at utsparingen ikke bryter gjennom kjernen på TR-1)
    b_kjerne=6.0,  # kjernebredde (snorbane). 7 -> 6 2026-11-05: TR-3 på 22,5 mm var 2 mm for lang
    hoyde_flens=3.0,  # flensens høyde over kjernen
    skraa=3.0,  # aksiell lengde av ytre flens' skrå innside (45°)
    t_ytre=0.5,  # rett del av ytre flens
    d_nav=20.0,
    l_nav=7.5,  # mutteren (6,6 over hjørnene) får ca. 0,45 mm vegg i hver ende
    slisse=1.0,
    e_skrue=6.0,  # skrueaksens avstand fra senter
    m3_hull=3.4,
    d_hode=6.2,  # forsenkning for skruehodet (linsehode Ø5,6)
    mutter=5.7,  # over flatene (M3-mutter 5,5)
    x_sete=4.25,  # hodesete og mutterflate, avstand fra slissen (klemlengde 8,5 for M3 x 12)
    krok_stilk=2.0,  # stilkens tykkelse i rotasjonsretningen
    krok_finger=2.0,  # hvor langt fingeren stikker ut bakover (mot rotasjonen)
    krok_gap=1.3,  # radiell klaring under fingeren, der løkka ligger (snor ca. 1 mm)
    krok_hoyde=2.3,  # krokens topp over kjerneflaten (flensen er 3,0)
    krok_klaring=1.0,  # aksiell klaring fra indre flens til krokens fot
    krok_bredde=4.0,  # aksiell bredde (de nederste ca. 2,5 mm er skrådd 45°), 1 mm klaring opp
)


def syl(r0, z0, z1, r1=None):
    return mf.Manifold.cylinder(z1 - z0, r0, r0 if r1 is None else r1, SEG).translate([0, 0, z0])


def boks(x0, x1, y0, y1, z0, z1):
    return mf.Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])


def lengder(p, d_kjerne):
    rk = d_kjerne / 2
    rf = rk + p["hoyde_flens"]
    z1 = p["t_indre"]
    z2 = z1 + p["b_kjerne"]
    z3 = z2 + p["skraa"]
    z4 = z3 + p["t_ytre"]
    z5 = z4 + p["l_nav"]
    return rk, rf, (z1, z2, z3, z4, z5)


def utsparing(p):
    """Plass til sveisekransen: sylinder Ø d_utsparing, deretter 45° kjegle inn til boringen."""
    ru, zr, rb = p["d_utsparing"] / 2, p["rett_utsparing"], p["boring"] / 2
    return syl(ru, -0.01, zr) + syl(ru, zr, zr + ru - rb + 0.3, rb - 0.3)


TH_KROK = math.pi  # krokens vinkelposisjon (på -x)


def krok_profil(p, r0=-1.0):
    """Krokprofilen i lokale koordinater (t, r): t langs rotasjonen under opptrekk, r utover fra
    kjerneflaten. Stilken er sentrert på t = 0, fingeren peker mot -t. r0 er fotens radius."""
    s, L, g, h = p["krok_stilk"] / 2, p["krok_finger"], p["krok_gap"], p["krok_hoyde"]
    return [(-s, r0), (s, r0), (s, h), (-s - L, h), (-s - L, g), (-s, g)]


def krok(p, rk, z1):
    c, si = math.cos(TH_KROK), math.sin(TH_KROK)
    er, et = (c, si), (-si, c)  # radiell og tangentiell (+theta) enhetsvektor
    pol = [((rk + r) * er[0] + t * et[0], (rk + r) * er[1] + t * et[1]) for t, r in krok_profil(p)]
    if sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pol, pol[1:] + pol[:1])) < 0:
        pol.reverse()
    zs = z1 + p["krok_klaring"]
    ze = zs + p["krok_bredde"]
    k = mf.Manifold.extrude(mf.CrossSection([pol]), ze - zs).translate([0, 0, zs])
    # 45° skrå underside: kroken vokser ut fra kjernen, så ingen del av den henger fritt
    return k ^ syl(rk, zs, ze, rk + (ze - zs))


def klemnav_kutt(p, rb, z4, z5):
    """Det som skal trekkes fra et Ø d_nav-nav på z4..z5: to slisser (±y) og to M3-skruer langs x
    med hodesete og mutterlomme, rotert 180° fra hverandre for balanse. Brukes også av pinjong.py."""
    s = p["slisse"] / 2
    zm = (z4 + z5) / 2
    rn = p["d_nav"] / 2 + 1
    feste = boks(-s, s, rb - 0.5, rn, z4, z5 + 1)
    hull = mf.Manifold.cylinder(2 * rn, p["m3_hull"] / 2, p["m3_hull"] / 2, 48).rotate([0, 90, 0])
    hull = hull.translate([-rn, p["e_skrue"], zm])
    hode = mf.Manifold.cylinder(rn, p["d_hode"] / 2, p["d_hode"] / 2, 48).rotate([0, -90, 0])
    hode = hode.translate([-p["x_sete"], p["e_skrue"], zm])  # fra hodesetet og utover mot -x
    r_hex = p["mutter"] / math.sqrt(3)  # omskrevet radius
    hex_ = mf.Manifold.extrude(mf.CrossSection.circle(r_hex, 6), rn)  # spiss i ±x i snittet
    hex_ = hex_.rotate([0, 90, 0])  # aksen langs +x; snittet ligger i yz med spiss mot ±z
    hex_ = hex_.translate([p["x_sete"], p["e_skrue"], zm])
    feste = feste + hull + hode + hex_
    return feste + feste.rotate([0, 0, 180])


def trommel(p, d_kjerne):
    rk, rf, (z1, z2, z3, z4, z5) = lengder(p, d_kjerne)
    rb = p["boring"] / 2

    # Rotasjonsdelen: indre flens, kjerne, skrå + rett ytre flens, nav
    m = syl(rf, 0, z1) + syl(rk, z1, z2) + syl(rk, z2, z3, rf) + syl(rf, z3, z4) + syl(p["d_nav"] / 2, z4, z5)
    m = m - syl(rb, -1, z5 + 1) - utsparing(p)

    m = m + krok(p, rk, z1)

    return m - klemnav_kutt(p, rb, z4, z5)


def skriv_3mf(m, sti, navn):
    mesh = m.to_mesh()
    v = np.asarray(mesh.vert_properties)[:, :3]
    f = np.asarray(mesh.tri_verts)
    verts = "".join(f'<vertex x="{a:.4f}" y="{b:.4f}" z="{c:.4f}"/>' for a, b, c in v)
    tris = "".join(f'<triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in f)
    modell = ('<?xml version="1.0" encoding="UTF-8"?>\n'
              '<model unit="millimeter" xml:lang="en-US" '
              'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">'
              f'<metadata name="Title">{navn}</metadata><resources>'
              f'<object id="1" type="model" name="{navn}"><mesh><vertices>{verts}</vertices>'
              f'<triangles>{tris}</triangles></mesh></object></resources>'
              '<build><item objectid="1"/></build></model>')
    typer = ('<?xml version="1.0" encoding="UTF-8"?>\n'
             '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
             '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>'
             '</Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Target="/3D/3dmodel.model" Id="rel0" '
            'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')
    with zipfile.ZipFile(sti, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", typer)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", modell)


def nokkeltall(p, d_kjerne, m):
    rk, rf, zs = lengder(p, d_kjerne)
    rb = p["boring"] / 2
    omkrets = 2 * math.pi * (rk + 0.5)  # snorsenter for ca. 1 mm snor
    return {
        "kjerne Ø": d_kjerne,
        "flens Ø": 2 * rf,
        "lengde": zs[-1],
        "y i bilen (fra-til)": (-3.0, -3.0 - zs[-1]),
        "vindinger for 20 cm": 200 / omkrets,
        "krok fra-til z": (zs[0] + p["krok_klaring"], zs[0] + p["krok_klaring"] + p["krok_bredde"]),
        "masse PLA 100 % [g]": m.volume() * RHO_PLA,
    }


def tegn(p, tromler, sti):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig, akser = plt.subplots(2, len(tromler), figsize=(4.2 * len(tromler), 8.4))
    for j, (navn, d, m) in enumerate(tromler):
        rk, rf, zs = lengder(p, d)
        # Snitt langs aksen i planet y = 0: går gjennom kroken (-x) og navet mellom skruene
        a = akser[0, j]
        snitt = m.rotate([90, 0, 0]).slice(0.0)  # (x, y, z) -> (x, -z, y)
        for poly in snitt.to_polygons():
            a.add_patch(Polygon(np.asarray(poly) * [1, -1], fc="#d4952f", ec="k", lw=0.6))
        a.set_title(f"{navn}: snitt langs aksen (y = 0)")
        a.set_xlabel("x [mm] (kroken på -x)"); a.set_ylabel("z [mm] (0 = mot svinghjulet)")
        a.set_xlim(-12.5, 12.5); a.set_ylim(-1, zs[-1] + 1); a.set_aspect("equal")
        # Sveisekransen Ø13 x 3 mm, tegnet som en konveks kvart-ellipse (verste tilfelle)
        rr = np.linspace(3.0, 6.5, 40)
        zk = 3.0 * np.sqrt(np.clip(1 - ((rr - 3.0) / 3.5) ** 2, 0, 1))
        for sgn in (-1, 1):
            a.plot(sgn * rr, zk, color="#2a5db5", lw=0.9, ls="--")
        a.text(0, 1.0, "sveis", ha="center", fontsize=7, color="#2a5db5")

        # Planriss sett fra navenden: kjernesnitt med kroken og navsnitt med skruene
        a = akser[1, j]
        for z, farge in (((zs[3] + zs[4]) / 2, "#d4952f"), (zs[0] + p["krok_klaring"] + p["krok_bredde"] - 0.4, "#f0d9a8")):
            for poly in m.slice(z).to_polygons():
                a.add_patch(Polygon(np.asarray(poly), fc=farge, ec="k", lw=0.6, alpha=0.85))
        th = np.linspace(math.pi * 0.25, math.pi * 0.75, 30)
        a.plot(1.15 * rf * np.cos(th), 1.15 * rf * np.sin(th), color="#b5322a", lw=1)
        a.annotate("", xy=(1.15 * rf * math.cos(th[-1]), 1.15 * rf * math.sin(th[-1])),
                   xytext=(1.15 * rf * math.cos(th[-3]), 1.15 * rf * math.sin(th[-3])),
                   arrowprops=dict(arrowstyle="->", color="#b5322a"))
        a.text(0, 1.22 * rf, "rotasjon under opptrekk", ha="center", fontsize=7, color="#b5322a")
        a.text(-rk - 2.6, -2.5, "krok", ha="right", fontsize=7)
        a.set_title(f"{navn}: sett fra navenden (+z)")
        a.set_xlim(-14, 14); a.set_ylim(-14, 16); a.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(sti, dpi=160)
    plt.close(fig)


def tegn_krok(p, d_kjerne, sti):
    """Detalj av kroken med løkka på, i to riss, for TR-1."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

    rk, rf, (z1, z2, z3, z4, z5) = lengder(p, d_kjerne)
    s, L, g, h = p["krok_stilk"], p["krok_finger"], p["krok_gap"], p["krok_hoyde"]
    zs, ze = z1 + p["krok_klaring"], z1 + p["krok_klaring"] + p["krok_bredde"]
    d_snor, rod, gul = 1.0, "#b5322a", "#d4952f"
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.6))

    # 1) Sett langs aksen: krokprofilen på kjerneflaten (rettet ut), løkka i krokens innerhjørne
    t = np.linspace(-6, 8, 50)
    a.fill_between(t, -2.5, 0, color="#f0d9a8", ec="k", lw=0.6)
    a.add_patch(Polygon(krok_profil(p, 0.0), fc=gul, ec="k", lw=0.8))
    s0 = -s / 2  # stilkens bakside (mot rotasjonen), der løkka ligger an
    a.axhline(p["hoyde_flens"], color="0.5", ls="--", lw=0.7)
    a.text(7.8, p["hoyde_flens"] + 0.1, "flenstopp", ha="right", va="bottom", fontsize=8, color="0.4")
    a.add_patch(plt.Circle((s0 - d_snor / 2 - 0.05, d_snor / 2), d_snor / 2, fc=rod, ec="k", lw=0.5))
    a.plot([s / 2, 7.5], [d_snor / 2] * 2, color=rod, lw=3, solid_capstyle="butt")
    a.annotate("", xy=(7.5, 1.3), xytext=(4.5, 1.3), arrowprops=dict(arrowstyle="->", color=rod))
    a.text(7.8, 1.6, "trekk under opplading\nog opptrekk", ha="right", va="bottom", fontsize=8, color=rod)
    a.annotate("løkka glir av her\nnår trommelen spinner videre", xy=(s0 - L - 0.3, h + 0.2),
               xytext=(-5.8, 3.9), fontsize=8, arrowprops=dict(arrowstyle="->", color="k"))
    a.text(s0 - L / 2, g + (h - g) / 2, "finger", ha="center", va="center", fontsize=7)
    a.text(0, h / 2 + 0.3, "stilk", ha="center", va="center", fontsize=7, rotation=90)
    a.text(-5.8, -1.4, "kjerne", fontsize=8)
    a.set_title("Sett langs aksen (kjerneflaten rettet ut)")
    a.set_xlabel("langs omkretsen [mm]  (+ = rotasjon under opptrekk)"); a.set_ylabel("radielt [mm]")
    a.set_xlim(-6, 8); a.set_ylim(-2.5, 5); a.set_aspect("equal")

    # 2) Sett radielt inn mot kjernen: løkka rundt stilken, fingeren over, flensene på sidene
    b.add_patch(Rectangle((-6, z1 - 1.5), 14, 1.5, fc=gul, ec="k", lw=0.6))
    b.add_patch(Rectangle((-6, z2), 14, 1.5, fc=gul, ec="k", lw=0.6, alpha=0.6))
    b.text(7.8, z1 - 0.75, "indre flens", ha="right", va="center", fontsize=8)
    b.text(7.8, z2 + 0.75, "ytre flens (skrå)", ha="right", va="center", fontsize=8)
    b.add_patch(Rectangle((-6, z1), 14, z2 - z1, fc="#f0d9a8", ec="none"))
    b.add_patch(Rectangle((s0 - L, zs), L + s, ze - zs, fc=gul, ec="k", lw=0.6, ls="--", alpha=0.45))
    b.add_patch(Rectangle((s0, zs), s, ze - zs, fc=gul, ec="k", lw=0.8))
    b.annotate("finger (over løkka)", xy=(s0 - L / 2, zs + 0.6), xytext=(-5.8, zs - 0.3), fontsize=7,
               va="center", arrowprops=dict(arrowstyle="-", color="k", lw=0.6))
    o = 0.15 + d_snor / 2
    b.add_patch(FancyBboxPatch((s0 - o, zs - o), s + 2 * o, ze - zs + 2 * o, boxstyle="round,pad=0,rounding_size=0.8",
                               fc="none", ec=rod, lw=3))
    b.plot([-s0 + o, 7.5], [(zs + ze) / 2] * 2, color=rod, lw=3)
    b.annotate("", xy=(7.5, (zs + ze) / 2 + 1.0), xytext=(5.0, (zs + ze) / 2 + 1.0),
               arrowprops=dict(arrowstyle="->", color=rod))
    b.text(-5.8, (zs + ze) / 2 + 0.8, "løkka ligger rundt\nstilkfoten, under\nfingeren", fontsize=8, va="center")
    b.set_title("Sett radielt inn mot kjernen")
    b.set_xlabel("langs omkretsen [mm]"); b.set_ylabel("z [mm] (0 = mot svinghjulet)")
    b.set_xlim(-6, 8); b.set_ylim(z1 - 1.5, z2 + 1.5); b.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(sti, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kjerner", type=float, nargs="+", default=[12.0, 14.0, 18.0], help="kjernediametre [mm]")
    ap.add_argument("--boring", type=float, default=P["boring"])
    a = ap.parse_args()
    p = dict(P, boring=a.boring)
    UT.mkdir(parents=True, exist_ok=True)
    tromler = []
    for i, d in enumerate(a.kjerner, 1):
        navn = f"TR-{i}"
        m = trommel(p, d)
        assert m.status() == mf.Error.NoError and m.genus() >= 0, navn
        skriv_3mf(m, UT / f"{navn}_kjerne{d:g}.3mf", f"Snortrommel {navn} kjerne Ø{d:g}")
        tromler.append((navn, d, m))
        print(navn, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in nokkeltall(p, d, m).items()})
    tegn(p, tromler, UT / "snortrommel_forhandsvisning.png")
    tegn_krok(p, tromler[0][1], UT / "snortrommel_krok.png")
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
