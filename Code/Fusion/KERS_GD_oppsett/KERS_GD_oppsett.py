"""Fusion-skript: bygger oppsettet for Generative Design av rammen til KERS-bilen.

Lager tre komponenter i et nytt design:
  GD Bevar     – preserve geometry: lagerhus (4 lukkede + 2 U-seter for svinghjulet) og
                 tverrbjelker foran og bak i full bredde (sikrer at rammen blir et rektangel)
  GD Hindring  – obstacle geometry: hjul, aksler, lagre, svinghjul, tannhjul, trommel,
                 snorbane og monteringsbaner (halvgjennomsiktige)
  GD Grenser   – obstacle geometry: 140 x 85 x 70 og bakkeklaring (skjul med øyet for å se inn)
  Lagerlokk    – to lokk som printes separat og skrus på U-setene (ikke med i studien)

Mål fra drivverk-tegningen: aksler Ø6 på z = 34, x = 34 / 70 / 106, 626-lagre (6 x 19 x 6)
med innerflate på y = ±26,5, Ø68 hjul, svinghjul Ø60 x 6, tannhjul z73 m0,8 (tupp Ø60).
x fra bakkant, y fra midtlinjen (+y = tannhjulssiden), z fra bakken. Mål i mm.

Kjør i Fusion: Utilities > Add-Ins > Scripts and Add-Ins > + (Script from device) > velg
mappen KERS_GD_oppsett > Run. Modellen bygges med Z opp.

Filen kan også importeres utenfor Fusion (geometri() bruker ikke adsk), slik at
gd_forhandsvisning.py kan tegne og kontrollere den samme geometrien.
"""

# ---------------------------------------------------------------- parametere [mm]
X_BAK, X_SV, X_FOR = 34.0, 70.0, 106.0
Z_A = 34.0
R_HJUL, B_HJUL, Y_HJUL = 34.0, 8.0, 33.5  # hjulradius, bredde, innerflate
D_LAGER, B_LAGER = 19.0, 6.0  # 626: ytre Ø og bredde (juster boringen etter testprint)
Y_INN = 26.5  # lagerets innerflate
R_HUS = 12.5  # ytre radius lagerhus (3 mm vegg)
LEPPE, D_LEPPE = 1.0, 17.0  # innvendig kant som holder lageret sideveis (treffer bare ytterringen)
L_SV_AKSEL = 65.0  # svinghjulsakselens lengde (oppgaveteksten). OBS: sjekk mot målt aksel!
LUFT = 1.0  # luft rundt roterende deler
M3_KJERNE = 2.8  # hull for M3 skrudd rett i PLA, 10 dypt (M3 x 16 gjennom 7 mm lokk)
M3_KLARING = 3.4
Z_MAKS, Y_MAKS, X_MAKS = 70.0, 42.5, 140.0
Z_BAKKE = 3.0  # minste bakkeklaring for rammen
# Innerrom: hindring som holder midten åpen, slik at rammen blir et rektangel sett ovenfra
# (sideskinner + tverrbjelke bak for x < INNERROM_X[0] og foran for x > INNERROM_X[1]).
INNERROM_X = (12.0, 112.0)
# Tverrbjelker i full bredde (y = ±32,5) som bevart geometri: x0, x1, z0, z1
BJELKE_FORAN = (120.0, 126.0, 41.0, 47.0)  # her holdes bilen under opptrekk
BJELKE_BAK = (3.0, 9.0, 4.0, 10.0)  # under tannhjulet og snorbanen

Y_LAGER = (Y_INN, Y_INN + B_LAGER)  # 26,5 – 32,5
Y_HUS = (Y_INN - LEPPE, Y_INN + B_LAGER)  # 25,5 – 32,5
R_L = D_LAGER / 2
# U-sete for svinghjulslageret
U_X = (X_SV - 16.0, X_SV + 16.0)  # 54 – 86
U_Z = (Z_A - R_L - 3.0, Z_A + 6.0)  # 21,5 – 40
X_SKRUER = (X_SV - 12.75, X_SV + 12.75)  # 57,25 og 82,75
Y_SKRUE = Y_INN + B_LAGER / 2  # 29,5
H_LOKK = 7.0


# ---------------------------------------------------------------- geometri (ren Python)
def boks(x0, x1, y0, y1, z0, z1):
    return ("boks", x0, x1, y0, y1, z0, z1)


def syl_y(x, z, r, y0, y1):
    return ("syl_y", x, z, r, y0, y1)


def syl_z(x, y, r, z0, z1):
    return ("syl_z", x, y, r, z0, z1)


def speil(prim):
    """Speiler en primitiv om y = 0."""
    t = prim[0]
    if t == "boks":
        _, x0, x1, y0, y1, z0, z1 = prim
        return boks(x0, x1, -y1, -y0, z0, z1)
    if t == "syl_y":
        _, x, z, r, y0, y1 = prim
        return syl_y(x, z, r, -y1, -y0)
    _, x, y, r, z0, z1 = prim
    return syl_z(x, -y, r, z0, z1)


def kropp(gruppe, navn, *ops):
    """ops: ("+", primitiv) eller ("-", primitiv). Første op må være "+"."""
    return dict(gruppe=gruppe, navn=navn, ops=list(ops))


def speil_kropp(k, navn):
    return dict(gruppe=k["gruppe"], navn=navn, ops=[(s, speil(p)) for s, p in k["ops"]])


def geometri():
    k = []
    y0, y1 = Y_LAGER
    h0, h1 = Y_HUS

    # ---- Bevar: lukkede lagerhus bak og foran (boring Ø19 + leppe Ø17 innerst)
    for x, navn in ((X_BAK, "Lagerhus bak"), (X_FOR, "Lagerhus foran")):
        hus = kropp("Bevar", f"{navn} +y",
                    ("+", syl_y(x, Z_A, R_HUS, h0, h1)),
                    ("-", syl_y(x, Z_A, R_L, y0, y1)),
                    ("-", syl_y(x, Z_A, D_LEPPE / 2, h0, y0)))
        k += [hus, speil_kropp(hus, f"{navn} -y")]

    # ---- Bevar: U-sete for svinghjulslageret, åpent oppover, med to skruehull for lokket
    u = kropp("Bevar", "U-sete svinghjul +y",
              ("+", boks(*U_X, h0, h1, *U_Z)),
              ("-", syl_y(X_SV, Z_A, R_L, y0, y1)),
              ("-", boks(X_SV - R_L, X_SV + R_L, y0, y1, Z_A, U_Z[1])),
              ("-", syl_y(X_SV, Z_A, D_LEPPE / 2, h0, y0)),
              ("-", boks(X_SV - D_LEPPE / 2, X_SV + D_LEPPE / 2, h0, y0, Z_A, U_Z[1])),
              *[("-", syl_z(xs, Y_SKRUE, M3_KJERNE / 2, U_Z[1] - 10, U_Z[1])) for xs in X_SKRUER])
    k += [u, speil_kropp(u, "U-sete svinghjul -y")]

    # ---- Bevar: tverrbjelker foran og bak i full bredde
    bx0, bx1, bz0, bz1 = BJELKE_FORAN
    k.append(kropp("Bevar", "Tverrbjelke foran", ("+", boks(bx0, bx1, -h1, h1, bz0, bz1))))
    bx0, bx1, bz0, bz1 = BJELKE_BAK
    k.append(kropp("Bevar", "Tverrbjelke bak", ("+", boks(bx0, bx1, -h1, h1, bz0, bz1))))

    # ---- Hindring: roterende deler med luft
    for x in (X_BAK, X_FOR):
        hjul = kropp("Hindring", f"Hjul x{x:.0f} +y", ("+", syl_y(x, Z_A, R_HJUL + LUFT, Y_HJUL - 0.5, Y_MAKS)))
        k += [hjul, speil_kropp(hjul, f"Hjul x{x:.0f} -y")]
        k.append(kropp("Hindring", f"Aksel x{x:.0f}", ("+", syl_y(x, Z_A, 3.5, -50, 50))))
    k.append(kropp("Hindring", "Svinghjulsaksel", ("+", syl_y(X_SV, Z_A, 3.5, -L_SV_AKSEL / 2, L_SV_AKSEL / 2))))
    for x in (X_BAK, X_SV, X_FOR):
        lag = kropp("Hindring", f"Lager x{x:.0f} +y", ("+", syl_y(x, Z_A, R_L, y0, y1)))
        k += [lag, speil_kropp(lag, f"Lager x{x:.0f} -y")]
    k.append(kropp("Hindring", "Svinghjul", ("+", syl_y(X_SV, Z_A, 30 + 2, -3 - 0.5, 3 + 0.5))))
    k.append(kropp("Hindring", "Tannhjul z73", ("+", syl_y(X_BAK, Z_A, 30 + LUFT, 4.5, 10.5))))
    k.append(kropp("Hindring", "Pinjong z17", ("+", syl_y(X_SV, Z_A, 7.6 + LUFT, 4.5, 10.5))))
    k.append(kropp("Hindring", "Nav tannhjul", ("+", syl_y(X_BAK, Z_A, 7 + LUFT, 9.5, 16.5))))
    k.append(kropp("Hindring", "Nav pinjong", ("+", syl_y(X_SV, Z_A, 7 + LUFT, 9.5, 16.5))))
    k.append(kropp("Hindring", "Snortrommel", ("+", syl_y(X_SV, Z_A, 11 + LUFT, -14, -3.5))))

    # ---- Hindring: snorbane bakover fra toppen av trommelen
    k.append(kropp("Hindring", "Snorbane", ("+", boks(-5, X_SV, -14, -3.5, 37, 55))))

    # ---- Hindring: monteringsbaner (svinghjulet løftes inn ovenfra, tannhjulet likeså)
    top = Z_MAKS + 5
    lb = kropp("Hindring", "Bane lager +y", ("+", boks(X_SV - R_L, X_SV + R_L, y0, y1, Z_A, top)))
    k += [lb, speil_kropp(lb, "Bane lager -y")]
    k.append(kropp("Hindring", "Bane aksel", ("+", boks(X_SV - 3.5, X_SV + 3.5, -y0, y0, Z_A, top))))
    k.append(kropp("Hindring", "Bane pinjong og trommel", ("+", boks(X_SV - 12, X_SV + 12, -14, 16.5, Z_A, top))))
    k.append(kropp("Hindring", "Bane svinghjul", ("+", boks(X_SV - 30, X_SV + 30, -3.5, 3.5, Z_A, top))))
    k.append(kropp("Hindring", "Bane tannhjul", ("+", boks(X_BAK - 31, X_BAK + 31, 4.5, 10.5, Z_A, top))))
    k.append(kropp("Hindring", "Bane nav", ("+", boks(X_BAK - 8, X_BAK + 8, 9.5, 16.5, Z_A, top))))

    # ---- Hindring: lagerlokk og skrutrekkerbane over U-setene
    lr = kropp("Hindring", "Lokkrom +y", ("+", boks(*U_X, h0, h1, U_Z[1], U_Z[1] + H_LOKK + 1)))
    k += [lr, speil_kropp(lr, "Lokkrom -y")]
    for i, xs in enumerate(X_SKRUER):
        sk = kropp("Hindring", f"Skrutrekker {i + 1} +y", ("+", syl_z(xs, Y_SKRUE, 3.5, U_Z[1] + H_LOKK + 1, top)))
        k += [sk, speil_kropp(sk, f"Skrutrekker {i + 1} -y")]

    # ---- Hindring: innerrommet mellom sideskinnene holdes åpent
    k.append(kropp("Hindring", "Innerrom", ("+", boks(*INNERROM_X, -h0, h0, Z_BAKKE, Z_MAKS + 5))))

    # ---- Hindring: grenser 140 x 85 x 70 og bakkeklaring
    k.append(kropp("Grenser", "Grense bakke", ("+", boks(-10, 150, -50, 50, -10, Z_BAKKE))))
    k.append(kropp("Grenser", "Grense tak", ("+", boks(-10, 150, -50, 50, Z_MAKS, Z_MAKS + 10))))
    k.append(kropp("Grenser", "Grense bak", ("+", boks(-10, 0, -50, 50, Z_BAKKE, Z_MAKS))))
    k.append(kropp("Grenser", "Grense foran", ("+", boks(X_MAKS, 150, -50, 50, Z_BAKKE, Z_MAKS))))
    k.append(kropp("Grenser", "Grense side +y", ("+", boks(0, X_MAKS, Y_MAKS, 50, Z_BAKKE, Z_MAKS))))
    k.append(kropp("Grenser", "Grense side -y", ("+", boks(0, X_MAKS, -50, -Y_MAKS, Z_BAKKE, Z_MAKS))))

    # ---- Lagerlokk (printes separat)
    lokk = kropp("Lagerlokk", "Lagerlokk +y",
                 ("+", boks(*U_X, h0, h1, U_Z[1], U_Z[1] + H_LOKK)),
                 ("-", syl_y(X_SV, Z_A, R_L, y0, y1)),
                 ("-", syl_y(X_SV, Z_A, D_LEPPE / 2, h0, y0)),
                 *[("-", syl_z(xs, Y_SKRUE, M3_KLARING / 2, U_Z[1], U_Z[1] + H_LOKK)) for xs in X_SKRUER])
    k += [lokk, speil_kropp(lokk, "Lagerlokk -y")]
    return k


# ---------------------------------------------------------------- Fusion
def _brep(tbm, prim, adsk):
    s = 0.1  # Fusion API bruker cm internt
    P = lambda x, y, z: adsk.core.Point3D.create(x * s, y * s, z * s)
    t = prim[0]
    if t == "boks":
        _, x0, x1, y0, y1, z0, z1 = prim
        obb = adsk.core.OrientedBoundingBox3D.create(
            P((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2),
            adsk.core.Vector3D.create(1, 0, 0), adsk.core.Vector3D.create(0, 1, 0),
            (x1 - x0) * s, (y1 - y0) * s, (z1 - z0) * s)
        return tbm.createBox(obb)
    if t == "syl_y":
        _, x, z, r, y0, y1 = prim
        return tbm.createCylinderOrCone(P(x, y0, z), r * s, P(x, y1, z), r * s)
    _, x, y, r, z0, z1 = prim
    return tbm.createCylinderOrCone(P(x, y, z0), r * s, P(x, y, z1), r * s)


def run(context):
    import traceback

    import adsk.core
    import adsk.fusion

    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        design = adsk.fusion.Design.cast(app.activeProduct)
        root = design.rootComponent
        tbm = adsk.fusion.TemporaryBRepManager.get()
        union = adsk.fusion.BooleanTypes.UnionBooleanType
        diff = adsk.fusion.BooleanTypes.DifferenceBooleanType

        kropper = geometri()
        grupper = {}
        for kr in kropper:
            grupper.setdefault(kr["gruppe"], []).append(kr)

        for gruppe, liste in grupper.items():
            occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
            comp = occ.component
            comp.name = "GD " + gruppe if gruppe != "Lagerlokk" else "Lagerlokk (printes separat)"
            gjennomsiktig = 0.3 if gruppe in ("Hindring", "Grenser") else 1.0
            base = comp.features.baseFeatures.add()
            base.startEdit()
            for kr in liste:
                (_, p0), *resten = kr["ops"]
                b = _brep(tbm, p0, adsk)
                for tegn, p in resten:
                    tbm.booleanOperation(b, _brep(tbm, p, adsk), union if tegn == "+" else diff)
                ny = comp.bRepBodies.add(b, base)
                ny.name = kr["navn"]
            base.finishEdit()
            if gjennomsiktig < 1.0:
                try:
                    occ.opacity = gjennomsiktig
                except Exception:
                    pass
            if gruppe == "Grenser":
                occ.isLightBulbOn = False  # skjult; slå på igjen før studien hvis Fusion krever det

        app.activeViewport.fit()
        # Ingen meldingsboks når alt går bra: på Mac kan den havne bak vinduet og låse menyene.
        try:
            app.log(f"KERS_GD_oppsett: ferdig, {len(kropper)} kropper. Preserve = GD Bevar, "
                    "Obstacle = GD Hindring + GD Grenser.")
        except Exception:
            pass
    except Exception:
        ui.messageBox("Feil:\n" + traceback.format_exc())
