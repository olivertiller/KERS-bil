"""Testringer for 626-lagrene i den 3D-printede rammen.

Hver ring er en kopi av lagerhusene i rammen: boring D (6 mm dyp), en 1 mm kant med Ø17 innerst
og ytre Ø25. Ringene merkes med hakk i ytterkanten: 1 hakk = første diameter, 2 hakk = andre osv.

    python3 Code/lagertestringer.py                       # Ø18,9–19,3
    python3 Code/lagertestringer.py --d 19.0 19.1 19.2    # egne diametre

Lager Bygg/Generativ_ramme/lagertestringer.3mf (mm) og én STL per ring.
Print med boringsåpningen ned mot platen og kanten opp, slik lagerhusene ligger når rammen
printes på siden. Bruk samme filament og profil som rammen.
"""

import argparse
import math
import zipfile
from pathlib import Path

import numpy as np

UT = Path(__file__).resolve().parent.parent / "Bygg" / "Generativ_ramme"
R_YTRE, B_LAGER, LEPPE, D_LEPPE = 12.5, 6.0, 1.0, 17.0
N = 360  # punkter rundt


def ring(d_boring, n_hakk):
    """Lukket trekantnett for én ring. Returnerer (punkter, trekanter)."""
    ri, rl = d_boring / 2, D_LEPPE / 2
    h = B_LAGER + LEPPE
    # Ytre radius med V-hakk (1 mm dype, 8° brede) med 18° mellomrom
    th = np.linspace(0, 2 * np.pi, N, endpoint=False)
    ro = np.full(N, R_YTRE)
    for k in range(n_hakk):
        senter = math.radians(90 + 18 * (k - (n_hakk - 1) / 2))
        dv = np.abs((th - senter + np.pi) % (2 * np.pi) - np.pi)
        ro = np.minimum(ro, R_YTRE - 1.0 * np.clip(1 - dv / math.radians(4), 0, 1))
    # Tverrsnitt (r, z) mot klokka sett i rz-planet: bunn, utside, topp, kant, boring
    prof = [("i", ri, 0.0), ("o", None, 0.0), ("o", None, h), ("l", rl, h), ("l", rl, B_LAGER), ("i", ri, B_LAGER)]
    P = len(prof)
    pts = []
    for j in range(N):
        c, s = math.cos(th[j]), math.sin(th[j])
        for typ, r, z in prof:
            rr = ro[j] if typ == "o" else r
            pts.append((rr * c, rr * s, z))
    tris = []
    for j in range(N):
        a, b = j * P, ((j + 1) % N) * P
        for p in range(P):
            q = (p + 1) % P
            # firkant mellom profilkant p→q på vinkel j og j+1, normal utover fra materialet
            tris += [(a + p, b + p, b + q), (a + p, b + q, a + q)]
    return np.array(pts), np.array(tris)


def sjekk(pts, tris):
    kanter = {}
    for t in tris:
        for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            kanter[e] = kanter.get(e, 0) + 1
    ok = all(kanter.get((b, a), 0) == 1 for (a, b), n in kanter.items() if n == 1) and all(n == 1 for n in kanter.values())
    v = sum(np.dot(pts[t[0]], np.cross(pts[t[1]], pts[t[2]])) for t in tris) / 6
    return ok, v


def skriv_stl(fil, pts, tris, navn):
    with open(fil, "w") as f:
        f.write(f"solid {navn}\n")
        for t in tris:
            a, b, c = pts[t]
            n = np.cross(b - a, c - a)
            n = n / (np.linalg.norm(n) or 1)
            f.write(f" facet normal {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n  outer loop\n")
            for p in (a, b, c):
                f.write(f"   vertex {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n")
            f.write("  endloop\n endfacet\n")
        f.write(f"endsolid {navn}\n")


def skriv_3mf(fil, objekter):
    obj_xml, items = [], []
    for i, (navn, pts, tris, dx) in enumerate(objekter, start=1):
        v = "".join(f'<vertex x="{p[0] + dx:.5f}" y="{p[1]:.5f}" z="{p[2]:.5f}"/>' for p in pts)
        t = "".join(f'<triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in tris)
        obj_xml.append(f'<object id="{i}" name="{navn}" type="model"><mesh><vertices>{v}</vertices>'
                       f'<triangles>{t}</triangles></mesh></object>')
        items.append(f'<item objectid="{i}"/>')
    model = ('<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xml:lang="en-US" '
             'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">'
             f'<resources>{"".join(obj_xml)}</resources><build>{"".join(items)}</build></model>')
    types = ('<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
             '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Target="/3D/3dmodel.model" Id="rel0" '
            'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')
    with zipfile.ZipFile(fil, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", types)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", model)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d", type=float, nargs="+", default=[18.9, 19.0, 19.1, 19.2, 19.3], help="boringsdiametre [mm]")
    a = ap.parse_args()
    UT.mkdir(parents=True, exist_ok=True)
    objekter = []
    for i, d in enumerate(a.d):
        pts, tris = ring(d, i + 1)
        ok, vol = sjekk(pts, tris)
        navn = f"Lagerring D{d:.2f} ({i + 1} hakk)"
        skriv_stl(UT / f"lagerring_D{d:.2f}.stl", pts, tris, navn)
        objekter.append((navn, pts, tris, 30.0 * i))
        print(f"{navn:<28} lukket={ok}  volum {vol:6.0f} mm³  ({vol * 1.24e-3:.1f} g PLA)")
    skriv_3mf(UT / "lagertestringer.3mf", objekter)
    print(f"Lagret {UT / 'lagertestringer.3mf'} og {len(objekter)} STL-filer")


if __name__ == "__main__":
    main()
