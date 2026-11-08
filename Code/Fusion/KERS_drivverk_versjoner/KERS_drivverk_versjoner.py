"""Fusion-skript: bygger de fire CAD-versjonene av drivverket (seksjon 5.1, tab:dv-samlet).

For hver versjon lages et nytt design med aksler (Ø3 sveisetråd eller printet Ø6 med 4 x 4-firkant i endene), lagre, svinghjul, tannhjulspar (med tenner) og hjul. Skriptet lagrer
  Bygg/Drivverk_versjoner/drivverk_v<n>.f3d
  Rapport/Images/Drivverk/drivverk_v<n>_side.png   (sett fra tannhjulssiden)
  Rapport/Images/Drivverk/drivverk_v<n>_skraatt.png
I sidebildet er hjul, lagre og hylser på tannhjulssiden skjult, i skråbildet bare hjulene.

Koordinater som i KERS_GD_oppsett: x fra bakkant, y fra midtlinjen (+y = tannhjulssiden), z fra bakken, mm.

Kjør: Utilities > Add-Ins > Scripts and Add-Ins > + > velg mappen KERS_drivverk_versjoner > Run.
Endre VERSJONER_SOM_BYGGES hvis bare noen skal bygges.

Filen kan importeres utenfor Fusion: python3 KERS_drivverk_versjoner.py skriver kontrollmålene.
"""

import math
import os

VERSJONER_SOM_BYGGES = (1, 2, 3, 4)

M = 0.8  # modul
Z_PINJONG = 17
R_SV = 30.0  # svinghjul Ø60 x 6
Y_LAGER = (26.5, 32.5)  # 626-lager 6 x 19 x 6
R_LAGER = 9.5
Y_HJUL = (33.5, 41.5)  # hjulbredde 8
L_SV_AKSEL = 65.0
Y_PINJONG = (5.0, 12.0)
Y_TANNHJUL = (5.0, 10.0)
Y_FIRKANT = 33.0  # firkantenden på hjulakslene starter her

# Versjonene i tab:dv-samlet. dx = avstand i x mellom akslene. I versjon 1 ligger svinghjulsakselen
# 4 mm høyere enn hjulakslene, så dx er satt slik at senteravstanden blir nøyaktig a = 34 mm.
VERSJONER = {
    1: dict(navn="Versjon 1: konsept 3", z_sv=36.0, z_hjul=32.0, r_hjul=32.0, d_aksel=3.0, z_tannhjul=68),
    2: dict(navn="Versjon 2: aksler i lik høyde", z_sv=34.0, z_hjul=34.0, r_hjul=34.0, d_aksel=3.0, z_tannhjul=68),
    3: dict(navn="Versjon 3: Ø6-aksler", z_sv=34.0, z_hjul=34.0, r_hjul=34.0, d_aksel=6.0, z_tannhjul=68),
    4: dict(navn="Versjon 4: a = 36 mm", z_sv=34.0, z_hjul=34.0, r_hjul=34.0, d_aksel=6.0, z_tannhjul=73),
}


def mal(v):
    """Avledede mål for én versjon."""
    p = dict(VERSJONER[v])
    a = M * (Z_PINJONG + p["z_tannhjul"]) / 2
    dz = p["z_sv"] - p["z_hjul"]
    dx = math.sqrt(a * a - dz * dz)
    p.update(a=a, dx=dx, x_bak=p["r_hjul"], x_sv=p["r_hjul"] + dx, x_for=p["r_hjul"] + 2 * dx)
    p["lengde"] = p["x_for"] + p["r_hjul"]
    p["klaring_aksel"] = a - R_SV - p["d_aksel"] / 2
    p["bakkeklaring"] = p["z_sv"] - R_SV
    p["i"] = p["z_tannhjul"] / Z_PINJONG
    return p


# ---------------------------------------------------------------- geometri (ren Python)
def syl_y(x, z, r, y0, y1):
    return ("syl_y", x, z, r, y0, y1)


def boks(x0, x1, y0, y1, z0, z1):
    return ("boks", x0, x1, y0, y1, z0, z1)


def tann(cx, cz, vinkel, r0, r1, tykkelse, y0, y1):
    """Radiell blokk fra radius r0 til r1 i retning vinkel (i xz-planet), tykkelse på tvers."""
    return ("tann", cx, cz, vinkel, r0, r1, tykkelse, y0, y1)


def tannhjul(cx, cz, z, vinkel0, y0, y1):
    """ops for et tannhjul med modul M og z tenner; første tann peker i retning vinkel0."""
    r_rot, r_del, r_tupp = M * (z - 2.5) / 2, M * z / 2, M * (z + 2) / 2
    ops = [("+", syl_y(cx, cz, r_rot + 0.05, y0, y1))]
    lag = ((r_rot, r_del - 0.4, 1.6 * M * 1.25), (r_del - 0.5, r_del + 0.4, 1.25 * 1.25 * M), (r_del + 0.3, r_tupp, 0.8 * 1.25 * M))
    for n in range(z):
        v = vinkel0 + 2 * math.pi * n / z
        ops += [("+", tann(cx, cz, v, a, b, t, y0, y1)) for a, b, t in lag]
    return ops


def kropp(navn, farge, *ops):
    """Kropper med navn som slutter på +y ligger på tannhjulssiden, nærmest kameraet i bildene."""
    return dict(navn=navn, farge=farge, ops=list(ops))


def geometri(v):
    p = mal(v)
    k = []
    ra = p["d_aksel"] / 2
    zs, zh = p["z_sv"], p["z_hjul"]
    xb, xs, xf = p["x_bak"], p["x_sv"], p["x_for"]

    # Samme underlag i alle versjonene, så tilpasset visning gir lik målestokk i bildene.
    k.append(kropp("Underlag", "underlag", ("+", boks(-15, 160, -55, 55, -1, 0))))

    # Hjulaksler: Ø3 sveisetråd (v1, v2) eller printet Ø6 med firkant 4 x 4 x 10 i hver ende (v3, v4),
    # firkanten fra |y| = 33 og ut (Code/hjul.py).
    for x, navn in ((xb, "Bakaksel (drivaksel)"), (xf, "Foraksel")):
        if p["d_aksel"] < 6:
            k.append(kropp(navn, "stal", ("+", syl_y(x, zh, ra, -Y_HJUL[1], Y_HJUL[1]))))
        else:
            k.append(kropp(navn, "printet", ("+", syl_y(x, zh, ra, -Y_FIRKANT, Y_FIRKANT)),
                           ("+", boks(x - 2, x + 2, Y_FIRKANT, Y_FIRKANT + 10, zh - 2, zh + 2)),
                           ("+", boks(x - 2, x + 2, -Y_FIRKANT - 10, -Y_FIRKANT, zh - 2, zh + 2))))
    k.append(kropp("Svinghjulsaksel", "stal", ("+", syl_y(xs, zs, 3.0, -L_SV_AKSEL / 2, L_SV_AKSEL / 2))))

    # 626-lagre, med hylse Ø6/Ø3 der akselen er Ø3
    for x, z, r_inn in ((xb, zh, ra), (xs, zs, 3.0), (xf, zh, ra)):
        for s in (1, -1):
            y0, y1 = sorted((s * Y_LAGER[0], s * Y_LAGER[1]))
            k.append(kropp(f"Lager x{x:.0f} {'+' if s > 0 else '-'}y", "lager",
                           ("+", syl_y(x, z, R_LAGER, y0, y1)), ("-", syl_y(x, z, 3.0, y0, y1))))
            if r_inn < 3.0:
                k.append(kropp(f"Hylse x{x:.0f} {'+' if s > 0 else '-'}y", "hylse",
                               ("+", syl_y(x, z, 3.0, y0, y1)), ("-", syl_y(x, z, r_inn, y0, y1))))

    # Svinghjul og tannhjulspar med tenner. Tannen er tilnærmet med tre radielle blokker som smalner
    # mot tuppen (rot, rundt delesirkelen, tupp). Tennene er fasert slik at paret står i inngrep.
    k.append(kropp("Svinghjul", "svinghjul", ("+", syl_y(xs, zs, R_SV, -3, 3)), ("-", syl_y(xs, zs, 3.0, -3, 3))))
    phi = math.atan2(zh - zs, xb - xs)  # retning fra pinjongen mot tannhjulet
    k.append(kropp(f"Pinjong z{Z_PINJONG}", "pinjong",
                   *tannhjul(xs, zs, Z_PINJONG, phi, *Y_PINJONG), ("-", syl_y(xs, zs, 3.0, *Y_PINJONG))))
    zt = p["z_tannhjul"]
    k.append(kropp(f"Tannhjul z{zt}", "tannhjul",
                   *tannhjul(xb, zh, zt, phi + math.pi + math.pi / zt, *Y_TANNHJUL), ("-", syl_y(xb, zh, ra, *Y_TANNHJUL))))

    # Hjul
    for x in (xb, xf):
        for s in (1, -1):
            y0, y1 = sorted((s * Y_HJUL[0], s * Y_HJUL[1]))
            k.append(kropp(f"Hjul x{x:.0f} {'+' if s > 0 else '-'}y", "hjul",
                           ("+", syl_y(x, zh, p["r_hjul"], y0, y1)),
                           ("-", syl_y(x, zh, ra, y0, y1) if p["d_aksel"] < 6 else boks(x - 2, x + 2, y0, y1, zh - 2, zh + 2)),
                           ))
    return k


FARGER = {  # RGB
    "stal": (150, 150, 155), "printet": (225, 200, 150), "lager": (90, 90, 95), "hylse": (200, 160, 60),
    "svinghjul": (70, 75, 85), "pinjong": (40, 110, 200), "tannhjul": (40, 110, 200),
    "hjul": (35, 35, 35), "underlag": (235, 235, 235),
}


# ---------------------------------------------------------------- Fusion
def _brep(tbm, prim, adsk):
    s = 0.1  # Fusion API bruker cm
    P = lambda x, y, z: adsk.core.Point3D.create(x * s, y * s, z * s)
    if prim[0] == "boks":
        _, x0, x1, y0, y1, z0, z1 = prim
        obb = adsk.core.OrientedBoundingBox3D.create(
            P((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2),
            adsk.core.Vector3D.create(1, 0, 0), adsk.core.Vector3D.create(0, 1, 0),
            (x1 - x0) * s, (y1 - y0) * s, (z1 - z0) * s)
        return tbm.createBox(obb)
    if prim[0] == "tann":
        _, cx, cz, v, r0, r1, t, y0, y1 = prim
        c, sn, rm = math.cos(v), math.sin(v), (r0 + r1) / 2
        obb = adsk.core.OrientedBoundingBox3D.create(
            P(cx + rm * c, (y0 + y1) / 2, cz + rm * sn),
            adsk.core.Vector3D.create(c, 0, sn), adsk.core.Vector3D.create(-sn, 0, c),
            (r1 - r0) * s, t * s, (y1 - y0) * s)
        return tbm.createBox(obb)
    _, x, z, r, y0, y1 = prim
    return tbm.createCylinderOrCone(P(x, y0, z), r * s, P(x, y1, z), r * s)


def _utseender(app, design, adsk):
    """Lager én farget kopi av et plastutseende per farge. Feiler stille (da blir alt grått)."""
    ut = {}
    try:
        lib = None
        for navn in ("Fusion Appearance Library", "Fusion 360 Appearance Library"):
            lib = app.materialLibraries.itemByName(navn)
            if lib:
                break
        if not lib:
            return ut
        grunn = None
        for i in range(lib.appearances.count):
            a = lib.appearances.item(i)
            if "Plastic" in a.name and "Matte" in a.name:
                grunn = a
                break
        if grunn is None:
            return ut
        for navn, (r, g, b) in FARGER.items():
            ap = design.appearances.itemByName("KERS " + navn) or design.appearances.addByCopy(grunn, "KERS " + navn)
            for j in range(ap.appearanceProperties.count):
                prop = ap.appearanceProperties.item(j)
                if prop.objectType == adsk.core.ColorProperty.classType():
                    try:
                        prop.value = adsk.core.Color.create(r, g, b, 255)
                    except Exception:
                        pass
            ut[navn] = ap
    except Exception:
        pass
    return ut


def _kamera(app, adsk, eye, target):
    s = 0.1
    vp = app.activeViewport
    cam = vp.camera
    cam.isSmoothTransition = False
    cam.cameraType = adsk.core.CameraTypes.OrthographicCameraType
    cam.target = adsk.core.Point3D.create(*(t * s for t in target))
    cam.eye = adsk.core.Point3D.create(*(e * s for e in eye))
    cam.upVector = adsk.core.Vector3D.create(0, 0, 1)
    cam.isFitView = True
    vp.camera = cam
    vp.refresh()
    adsk.doEvents()
    return vp


def bygg(app, adsk, v, repo):
    import adsk.fusion

    p = mal(v)
    app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    design = adsk.fusion.Design.cast(app.activeProduct)
    root = design.rootComponent
    tbm = adsk.fusion.TemporaryBRepManager.get()
    union = adsk.fusion.BooleanTypes.UnionBooleanType
    diff = adsk.fusion.BooleanTypes.DifferenceBooleanType
    farger = _utseender(app, design, adsk)

    base = root.features.baseFeatures.add()
    base.startEdit()
    for kr in geometri(v):
        (_, p0), *resten = kr["ops"]
        b = _brep(tbm, p0, adsk)
        for tegn, prim in resten:
            tbm.booleanOperation(b, _brep(tbm, prim, adsk), union if tegn == "+" else diff)
        ny = root.bRepBodies.add(b, base)
        ny.name = kr["navn"]
    base.finishEdit()

    # Kroppene får nye referanser etter finishEdit; sett navn/farge via navnet.
    for kr in geometri(v):
        body = root.bRepBodies.itemByName(kr["navn"])
        if body and kr["farge"] in farger:
            try:
                body.appearance = farger[kr["farge"]]
            except Exception:
                pass

    bildemappe = os.path.join(repo, "Rapport", "Images", "Drivverk")
    f3dmappe = os.path.join(repo, "Bygg", "Drivverk_versjoner")
    os.makedirs(bildemappe, exist_ok=True)
    os.makedirs(f3dmappe, exist_ok=True)

    midt = (72.5, 0.0, 30.0)  # likt i alle versjonene
    def vis(prefikser, synlig):
        """Slår av/på kroppene på tannhjulssiden (+y) som starter med et av prefiksene."""
        for kr in geometri(v):
            if kr["navn"].endswith("+y") and kr["navn"].startswith(prefikser):
                body = root.bRepBodies.itemByName(kr["navn"])
                if body:
                    body.isLightBulbOn = synlig
        adsk.doEvents()

    # Sidebilde: hjul, lagre og hylser nærmest kameraet skjult, så pinjongen og klaringen mellom
    # svinghjulet og forakselen synes.
    naer = ("Hjul", "Lager", "Hylse")
    vis(naer, False)
    vp = _kamera(app, adsk, (midt[0], 500.0, midt[2]), midt)
    vp.saveAsImageFile(os.path.join(bildemappe, f"drivverk_v{v}_side.png"), 1800, 1000)
    # Skråbilde forfra, ovenfra og fra tannhjulssiden, bare de nære hjulene skjult
    vis(("Lager", "Hylse"), True)
    vp = _kamera(app, adsk, (midt[0] + 260.0, 320.0, midt[2] + 220.0), midt)
    vp.saveAsImageFile(os.path.join(bildemappe, f"drivverk_v{v}_skraatt.png"), 1800, 1200)
    vis(naer, True)  # f3d-filen lagres med alt synlig

    em = design.exportManager
    em.execute(em.createFusionArchiveExportOptions(os.path.join(f3dmappe, f"drivverk_v{v}.f3d")))


def run(context):
    import traceback

    import adsk.core

    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        repo = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
        for v in VERSJONER_SOM_BYGGES:
            bygg(app, adsk, v, repo)
        try:
            app.log(f"KERS_drivverk_versjoner: ferdig med versjon {VERSJONER_SOM_BYGGES}. "
                    f"Filer i {os.path.join(repo, 'Bygg', 'Drivverk_versjoner')}")
        except Exception:
            pass
    except Exception:
        ui.messageBox("Feil:\n" + traceback.format_exc())


if __name__ == "__main__":
    print(f"{'':4}{'i':>6}{'a':>7}{'dx':>7}{'lengde':>8}{'klaring':>9}{'bakke':>7}")
    for v in VERSJONER:
        p = mal(v)
        print(f"v{v:<3}{p['i']:6.2f}{p['a']:7.1f}{p['dx']:7.2f}{p['lengde']:8.1f}"
              f"{p['klaring_aksel']:9.1f}{p['bakkeklaring']:7.1f}  {len(geometri(v))} kropper")
