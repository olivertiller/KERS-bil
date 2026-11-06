"""Drivhjul til KERS-bilen som 3MF-filer (mm): felg i PLA, tynt dekk i TPU, og et bart hjul.

    python3 Code/hjul.py                  # standardmål
    python3 Code/hjul.py --firkant 4.3    # firkanthull (4,3 valgt etter firkanttesten)
    python3 Code/hjul.py --diametre 68    # bare Ø68

Krever manifold3d (som snortrommel.py). Skriver til Kers-bil/Hjul/:
  hjul_felg_TPU.3mf   felg Ø65 for TPU-dekket (4 stk)
  dekk_TPU.3mf        dekk som strekkes over felgen, ferdig Ø68 (4 stk)
  hjul_bart_D68.3mf   hjul uten dekk, Ø68 (DH-1 i testserien)
  *_D66.3mf           samme design 2 mm mindre (felg Ø63, ferdig Ø66), som tillegg
  firkanttest.3mf     lite stykke med firkanthull 4,0–4,3 for å finne passningen
  hjul_forhandsvisning.png

Mål fra drivverk-tegningen: ytre Ø68 med dekk, hjulet på y ±33,5 til ±41,5 (8 bredt), 1 mm luft
til veggen (ytterflate y ±32,5). Akselendene freses til 4 x 4 mm firkant (diagonal 5,66 < Ø6).

Hjulet printes liggende med utsiden ned. Koordinater for delen: z = 0 er utsiden av hjulet,
z = 8 innsiden mot veggen. I bilen: y = ±(41,5 - z).
  * Nav Ø12 med firkanthull, ekstra hjørneavlastning så skarpe fresehjørner går inn.
  * Krage Ø7,8 x 1 på innsiden, ligger an mot innerringen i lageret (aldri ytterringen, Ø8 maks)
    og holder akselen på plass sideveis. Kragen og 0,5 mm av navet har rundt hull Ø6,4, så
    den runde delen av akselen like ved lageret får plass. Firkanten freses fra y ±33,0 og ut.
  * Seks eiker i en 2 mm skive på utsiden holder treghetsmomentet nede.
"""

import argparse
import math
from pathlib import Path

import manifold3d as mf
import numpy as np

from snortrommel import RHO_PLA, boks, skriv_3mf, syl

UT = Path(__file__).resolve().parent.parent / "Kers-bil" / "Hjul"
RHO_TPU = 1.21e-3  # g/mm^3
DIAMETRE = [68.0, 66.0]  # ferdige diametre som lages; Ø66 er et tillegg, 2 mm mindre

P = dict(
    d_ytre=68.0,  # ferdig hjul med dekk
    bredde=8.0,
    t_dekk=1.5,  # TPU-dekkets tykkelse
    b_dekk=6.0,  # dekkbredde (sporet mellom skuldrene)
    strekk=0.6,  # dekkets indre Ø er så mye mindre enn felgsetet
    skulder=0.8,  # skulderhøyde over setet
    t_felg=2.0,  # felgring, radielt
    t_skive=2.0,  # eikeskive
    n_eiker=6,
    b_eike=4.0,
    d_nav=12.0,
    firkant=4.3,  # tegnet firkanthull, valgt etter firkanttesten 2026-11-05 (akselen 4,00)
    d_krage=7.8,
    l_krage=1.0,
    d_rund=6.4,  # rundt hull i kragen for akselens runde del
    l_rund=1.5,  # fra innsiden (krage + 0,5 inn i navet)
    fas=0.5,
)


def firkanthull(s, z0, z1):
    """Firkanthull med hjørneavlastning (små sirkler i hjørnene)."""
    h = mf.Manifold.cube([s, s, z1 - z0], center=True).translate([0, 0, (z0 + z1) / 2])
    r = 0.45
    for sx in (-1, 1):
        for sy in (-1, 1):
            c = s / 2 - r * 0.35
            h = h + mf.Manifold.cylinder(z1 - z0, r, r, 24).translate([sx * c, sy * c, z0])
    return h


def eikeskive(p, r_inn, r_ut):
    t = p["t_skive"]
    skive = syl(r_ut, 0, t) - syl(r_inn, -1, t + 1)
    eiker = mf.Manifold()
    for i in range(p["n_eiker"]):
        e = boks(0, r_ut + 0.5, -p["b_eike"] / 2, p["b_eike"] / 2, 0, t)
        eiker = eiker + e.rotate([0, 0, 360 * i / p["n_eiker"]])
    return (syl(r_ut, 0, t) - skive) + (skive ^ eiker)  # hel ytterkant fjernes, eiker igjen


def hjul(p, med_dekk=True):
    b = p["bredde"]
    r_sete = p["d_ytre"] / 2 - p["t_dekk"] if med_dekk else p["d_ytre"] / 2
    r_felg_inn = r_sete - p["t_felg"] - (0 if med_dekk else 1.5)
    fas = p["fas"]

    if med_dekk:
        kant = (b - p["b_dekk"]) / 2
        rs = r_sete + p["skulder"]
        felg = (syl(r_sete, 0, b) + syl(rs, 0, kant) + syl(rs, b - kant, b))
    else:  # bart hjul med fas på begge kanter
        felg = syl(r_sete - fas, 0, fas, r_sete) + syl(r_sete, fas, b - fas) + syl(r_sete, b - fas, b, r_sete - fas)
    felg = felg - syl(r_felg_inn, -1, b + 1)

    rn = p["d_nav"] / 2
    m = felg + eikeskive(p, rn - 0.5, r_felg_inn + 0.5) + syl(rn, 0, b) + syl(p["d_krage"] / 2, b, b + p["l_krage"])
    z_topp = b + p["l_krage"]
    m = m - firkanthull(p["firkant"], -1, z_topp + 1)
    m = m - syl(p["d_rund"] / 2, z_topp - p["l_rund"], z_topp + 1)
    return m


def dekk(p):
    r_sete = p["d_ytre"] / 2 - p["t_dekk"]
    ri = r_sete - p["strekk"] / 2
    ro = ri + p["t_dekk"]
    b = p["b_dekk"] - 0.2
    f = 0.4
    m = syl(ro - f, 0, f, ro) + syl(ro, f, b - f) + syl(ro, b - f, b, ro - f)
    return m - syl(ri, -1, b + 1)


def firkanttest(p):
    verdier = [4.0, 4.1, 4.2, 4.3]
    m = boks(0, 12 * len(verdier) + 2, 0, 14, 0, 5)
    for i, v in enumerate(verdier):
        m = m - firkanthull(v, -1, 6).translate([7 + 12 * i, 8, 0])
        # Merking: i streker ved siden av hullet (1 = 4,0 ... 4 = 4,3)
        for k in range(i + 1):
            m = m - boks(3 + 12 * i + 1.6 * k, 3.8 + 12 * i + 1.6 * k, 1, 3, 4.4, 6)
    return m


def treghet(m, rho):
    """J om z-aksen for et lukket mesh, i kg m^2 (rho i g/mm^3)."""
    mesh = m.to_mesh()
    v = np.asarray(mesh.vert_properties)[:, :3]
    f = np.asarray(mesh.tri_verts)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    vol = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0  # tetraedre mot origo
    j = 0.0
    for k in (0, 1):  # x og y
        s = a[:, k] + b[:, k] + c[:, k]
        q = a[:, k] ** 2 + b[:, k] ** 2 + c[:, k] ** 2
        j += np.sum(vol / 20.0 * (q + s * s))
    return j * rho * 1e-3 * 1e-6  # g mm^2 -> kg m^2


def tegn(p, deler, sti):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig, ax = plt.subplots(1, 3, figsize=(14, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    farger = {"felg": "#d9d9d9", "dekk": "#2b2b2b", "bart": "#bfbfbf"}
    # Snitt langs aksen (planet y = 0), felg + dekk montert og bart hjul
    for a, navn, sett in ((ax[1], "Felg + TPU-dekk", ("felg", "dekk")), (ax[2], "Bart hjul (DH-1)", ("bart",))):
        for k in sett:
            m = deler[k]
            if k == "dekk":  # vis dekket montert midt på felgen (forenklet, ustrukket)
                m = m.translate([0, 0, (p["bredde"] - p["b_dekk"] + 0.2) / 2])
            for poly in m.rotate([90, 0, 0]).slice(0.0).to_polygons():
                pts = np.asarray(poly) * [1, -1]
                a.add_patch(Polygon(pts[:, [1, 0]], fc=farger[k], ec="k", lw=0.5))
        a.axvline(p["bredde"], color="#b4541a", ls=":", lw=0.8)
        a.text(p["bredde"] + 0.3, 2, "innside\n(mot vegg)", fontsize=7, color="#b4541a")
        a.set_xlim(-1, 11); a.set_ylim(-35, 35); a.set_aspect("equal")
        a.set_xlabel("z [mm] (0 = utsiden)"); a.set_title(f"{navn}: snitt langs aksen")
    # Sett fra utsiden. Snittene har hull, så de tegnes som én sammensatt sti (hull blir åpne).
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch

    def flate(a, snitt, **kw):
        verts, koder = [], []
        for poly in snitt.to_polygons():
            pts = np.asarray(poly)
            verts += list(pts) + [pts[0]]
            koder += [MPath.MOVETO] + [MPath.LINETO] * (len(pts) - 1) + [MPath.CLOSEPOLY]
        a.add_patch(PathPatch(MPath(verts, koder), **kw))

    a = ax[0]
    flate(a, deler["felg"].slice(1.0), fc=farger["felg"], ec="k", lw=0.5)
    flate(a, deler["dekk"].slice(2.0), fc=farger["dekk"], ec="k", lw=0.5)
    a.set_xlim(-36, 36); a.set_ylim(-36, 36); a.set_aspect("equal")
    a.set_title("Felg med dekk, sett fra utsiden (snitt z = 1)")
    fig.tight_layout()
    fig.savefig(sti, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--firkant", type=float, default=P["firkant"], help="tegnet firkanthull [mm]")
    ap.add_argument("--strekk", type=float, default=P["strekk"], help="dekkets indre Ø mindre enn setet [mm]")
    ap.add_argument("--diametre", type=float, nargs="+", default=DIAMETRE, help="ferdige hjuldiametre [mm]")
    a = ap.parse_args()
    UT.mkdir(parents=True, exist_ok=True)

    test = firkanttest(dict(P, firkant=a.firkant))
    skriv_3mf(test, UT / "firkanttest.3mf", "firkanttest")
    for d in a.diametre:
        p = dict(P, firkant=a.firkant, strekk=a.strekk, d_ytre=d)
        tillegg = "" if d == P["d_ytre"] else f"_D{d:g}"  # Ø68 beholder de gamle filnavnene
        deler = {"felg": hjul(p, True), "dekk": dekk(p), "bart": hjul(p, False)}
        filer = {"felg": f"hjul_felg_TPU{tillegg}", "dekk": f"dekk_TPU{tillegg}", "bart": f"hjul_bart_D{d:g}"}
        for k, m in deler.items():
            assert m.status() == mf.Error.NoError, k
            skriv_3mf(m, UT / f"{filer[k]}.3mf", filer[k])

        j_felg, j_dekk, j_bart = treghet(deler["felg"], RHO_PLA), treghet(deler["dekk"], RHO_TPU), treghet(deler["bart"], RHO_PLA)
        print(f"--- Ø{d:g} ---")
        print(f"Felg:  {deler['felg'].volume() * RHO_PLA:5.1f} g   Dekk: {deler['dekk'].volume() * RHO_TPU:4.1f} g   "
              f"Bart hjul: {deler['bart'].volume() * RHO_PLA:5.1f} g   (100 % fyll)")
        print(f"J per hjul: felg + dekk {j_felg + j_dekk:.2e} kg m^2, bart {j_bart:.2e} kg m^2")
        print(f"J_hjul for fire hjul (til kers_sim.py): felg + dekk {4 * (j_felg + j_dekk):.2e}, bart {4 * j_bart:.2e}")
        r_sete = p["d_ytre"] / 2 - p["t_dekk"]
        print(f"Felgsete Ø{2 * r_sete:.1f}, dekk innvendig Ø{2 * r_sete - p['strekk']:.1f} "
              f"({p['strekk'] / (2 * r_sete) * 100:.1f} % strekk), ferdig Ø{p['d_ytre']:g}")
        tegn(p, deler, UT / f"hjul_forhandsvisning{tillegg}.png")
    print("Skrevet til", UT)


if __name__ == "__main__":
    main()
