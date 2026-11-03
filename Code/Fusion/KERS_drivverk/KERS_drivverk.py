"""Fusion-skript: setter inn drivverket til KERS-bilen i det åpne designet, for kollisjonskontroll.

Delene har nominelle mål (ingen luft), i samme koordinater som KERS_GD_oppsett og rammen fra
Create Design From Outcome: x fra bakkant, y fra midtlinjen (+y = tannhjulssiden), z fra bakken.

Kjør: åpne rammen, Utilities > Add-Ins > Scripts and Add-Ins > + > velg mappen KERS_drivverk > Run.
Deretter Inspect > Interference: velg rammen og komponenten «Drivverk (kontroll)», Compute.
Include Coincident Faces skal være AV (lagrene ligger nøyaktig i boringene).
"""

import importlib.util
import os

# Gjenbruk parametere og primitiver fra oppsettskriptet, så målene alltid er like.
_her = os.path.dirname(os.path.abspath(__file__))
_gd_fil = os.path.join(_her, "..", "KERS_GD_oppsett", "KERS_GD_oppsett.py")
if not os.path.exists(_gd_fil):
    raise FileNotFoundError("Finner ikke KERS_GD_oppsett.py ved siden av denne mappen: " + _gd_fil)
_spec = importlib.util.spec_from_file_location("kers_gd", _gd_fil)
gd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gd)

boks, syl_y, kropp, speil_kropp = gd.boks, gd.syl_y, gd.kropp, gd.speil_kropp
R_AKSEL = 3.0


def geometri():
    k = []
    za, (y0, y1), (h0, h1) = gd.Z_A, gd.Y_LAGER, gd.Y_HUS
    yh0, yh1 = gd.Y_HJUL, gd.Y_HJUL + gd.B_HJUL  # 33,5 – 41,5
    G = "Drivverk"

    # Aksler
    for x, navn in ((gd.X_BAK, "Bakaksel"), (gd.X_FOR, "Foraksel")):
        k.append(kropp(G, navn, ("+", syl_y(x, za, R_AKSEL, -yh1, yh1))))
    k.append(kropp(G, "Svinghjulsaksel", ("+", syl_y(gd.X_SV, za, R_AKSEL, -gd.L_SV_AKSEL / 2, gd.L_SV_AKSEL / 2))))

    # 626-lagre (6 x 19 x 6)
    for x in (gd.X_BAK, gd.X_SV, gd.X_FOR):
        lag = kropp(G, f"Lager x{x:.0f} +y", ("+", syl_y(x, za, gd.R_L, y0, y1)), ("-", syl_y(x, za, R_AKSEL, y0, y1)))
        k += [lag, speil_kropp(lag, f"Lager x{x:.0f} -y")]

    # Svinghjul og tannhjulspar (tuppdiametre), nav og snortrommel
    k.append(kropp(G, "Svinghjul", ("+", syl_y(gd.X_SV, za, 30, -3, 3)), ("-", syl_y(gd.X_SV, za, R_AKSEL, -3, 3))))
    k.append(kropp(G, "Tannhjul z73", ("+", syl_y(gd.X_BAK, za, 30, 5, 10)), ("-", syl_y(gd.X_BAK, za, R_AKSEL, 5, 10))))
    k.append(kropp(G, "Nav tannhjul", ("+", syl_y(gd.X_BAK, za, 7, 10, 16)), ("-", syl_y(gd.X_BAK, za, R_AKSEL, 10, 16))))
    k.append(kropp(G, "Pinjong z17", ("+", syl_y(gd.X_SV, za, 7.6, 5, 10)), ("-", syl_y(gd.X_SV, za, R_AKSEL, 5, 10))))
    k.append(kropp(G, "Nav pinjong", ("+", syl_y(gd.X_SV, za, 7, 10, 16)), ("-", syl_y(gd.X_SV, za, R_AKSEL, 10, 16))))
    k.append(kropp(G, "Snortrommel",
                   ("+", syl_y(gd.X_SV, za, 11, -5, -4)),  # indre flens
                   ("+", syl_y(gd.X_SV, za, 7, -12, -5)),  # kjerne
                   ("+", syl_y(gd.X_SV, za, 11, -13.5, -12)),  # ytre flens
                   ("-", syl_y(gd.X_SV, za, R_AKSEL, -13.5, -4))))

    # Hjul Ø68 x 8
    for x in (gd.X_BAK, gd.X_FOR):
        hj = kropp(G, f"Hjul x{x:.0f} +y", ("+", syl_y(x, za, gd.R_HJUL, yh0, yh1)), ("-", syl_y(x, za, R_AKSEL, yh0, yh1)))
        k += [hj, speil_kropp(hj, f"Hjul x{x:.0f} -y")]

    # Lagerlokkene fra oppsettet
    k += [dict(lk, gruppe=G) for lk in gd.geometri() if lk["gruppe"] == "Lagerlokk"]
    return k


def run(context):
    import traceback

    import adsk.core
    import adsk.fusion

    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        if not design:
            ui.messageBox("Åpne designet med rammen først (Design-arbeidsområdet).")
            return
        root = design.rootComponent
        tbm = adsk.fusion.TemporaryBRepManager.get()
        union = adsk.fusion.BooleanTypes.UnionBooleanType
        diff = adsk.fusion.BooleanTypes.DifferenceBooleanType

        occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
        comp = occ.component
        comp.name = "Drivverk (kontroll)"
        direkte = design.designType == adsk.fusion.DesignTypes.DirectDesignType
        base = None
        if not direkte:
            base = comp.features.baseFeatures.add()
            base.startEdit()
        for kr in geometri():
            (_, p0), *resten = kr["ops"]
            b = gd._brep(tbm, p0, adsk)
            for tegn, p in resten:
                tbm.booleanOperation(b, gd._brep(tbm, p, adsk), union if tegn == "+" else diff)
            ny = comp.bRepBodies.add(b, base) if base else comp.bRepBodies.add(b)
            ny.name = kr["navn"]
        if base:
            base.finishEdit()
        app.activeViewport.fit()
    except Exception:
        ui.messageBox("Feil:\n" + traceback.format_exc())
