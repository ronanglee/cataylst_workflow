"""Code taken and refactored from Tiaporn"""

import os
from pathlib import Path  # type: ignore

import numpy as np  # type: ignore
from activity_selectivity.merge_databases import merge  # type: ignore
from ase.db import connect  # type: ignore
from utils import read_config, retreive_sql_correction, run_logger  # type: ignore

# constant
kb = 8.617e-5  # Boltzmann constant in eV/K
h = 4.136e-15  # Planck constant in eV*s
T = 298.15  # Temperature in K
k = kb * T * np.log(10)
conc = 10**-6  # concentration of metal ion

# gas molecules at T = 298.15 K, p = 1 bar
# eV/molecule (EDFT + ZPE + H - TS + BEEF correction)
g_h2 = -7.128  # eV/molecule (Free energy + BEEF correction)
g_h2o = -12.869  # eV/molecule (Free energy + BEEF correction)
database_dir = Path(__file__).parent.parent.parent / "runs" / "databases"

ph = 0
u = 0.7

config = read_config()


def acid_stability(
    ph: float, u: float, m: str, g_mn4: float, h_master: dict, g_mn4_ooh: float
) -> tuple[float, float]:
    """Calculate the relative stability of the intermediates.

    Args:
        ph: The pH of the solution.
        u: The applied potential.
        m: The metal.
        g_mn4: The energy of the metal ion.
        h_master: The energies of the hydrogen intermediates.
        g_mn4_ooh: The energy of the metal ion in the OOH intermediate.

    Returns:
        b_rel: The relative stability of the bare surface.
        ooh_rel: The relative stability of the OOH intermediate.
    """
    db_metals = connect(os.path.join(database_dir, "e_m.db"))
    g_m = db_metals.get(metal=m).energy / len(db_metals.get_atoms(metal=m))  # ev/atom
    dg = [
        9999
    ] * 42  # initialise more than needed just in case an element has a lot of rxns
    g_ch0 = list(h_master.values())[0]
    g_ch1 = list(h_master.values())[1]
    g_ch2 = list(h_master.values())[2]
    g_ch3 = list(h_master.values())[3]
    g_ch4 = list(h_master.values())[4]
    # all these rxns can be found in the publised paper "Effects of electrolyte anion
    # adsorption on the activity and stability of single atom electrocatalysts: Patniboon, Hansen"
    # "DOI: 10.1063/5.0125654"
    if m == "Pd":
        g_m2 = g_m + 2 * (0.951 + (1 / 2) * 0.0592 * np.log10(conc))
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
    elif m == "Pt":
        g_m2 = g_m + 2 * (1.18 + (1 / 2) * 0.0592 * np.log10(conc))
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
    elif m == "Cr":
        g_m2 = g_m + 2 * (-0.913 + (1 / 2) * 0.0592 * np.log10(conc))
        g_m3 = g_m + 3 * (-0.744 + (1 / 3) * 0.0592 * np.log10(conc))
        g_moh2 = g_m3 + g_h2o - 0.5 * g_h2 + 3.81 * k
        g_hmo4_1 = g_m3 + 4 * g_h2o + 3 * (1.350) - (7 / 2) * g_h2
        g_mo4_2 = g_m3 + 4 * g_h2o + 3 * (1.477) - 4 * g_h2
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
        # 2. dissolution to m+3
        dg[6] = g_m3 + g_ch0 - g_mn4 - 3 * u
        dg[7] = g_m3 + g_ch1 - g_mn4 - 2 * u - 0.5 * g_h2 + 1 * k * ph
        dg[8] = g_m3 + g_ch2 - g_mn4 - 1 * u - 1.0 * g_h2 + 2 * k * ph
        dg[9] = g_m3 + g_ch3 - g_mn4 - 0 * u - 1.5 * g_h2 + 3 * k * ph
        dg[10] = g_m3 + g_ch4 - g_mn4 + 1 * u - 2.0 * g_h2 + 4 * k * ph
        # 3. dissolution to g_moh2
        dg[11] = (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph) + g_ch0 - g_mn4 - 3 * u
        dg[12] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[13] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[14] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[15] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 5. dissolution to g_hmo4_1
        dg[16] = (
            (g_hmo4_1 + (7 / 2) * g_h2 - 7 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 3 * u
        )
        dg[17] = (
            (g_hmo4_1 + (7 / 2) * g_h2 - 7 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[18] = (
            (g_hmo4_1 + (7 / 2) * g_h2 - 7 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[19] = (
            (g_hmo4_1 + (7 / 2) * g_h2 - 7 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[20] = (
            (g_hmo4_1 + (7 / 2) * g_h2 - 7 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 6. dissolution to g_mo4_2
        dg[21] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 3 * u
        )
        dg[22] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[23] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[24] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[25] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
    elif m == "Mn":
        g_m2 = g_m + 2 * (-1.185 + (1 / 2) * 0.0592 * np.log10(conc))
        g_m3 = g_m2 + 1.5415
        g_mo4_2 = g_m2 + 4 * g_h2o - 4 * g_h2 + 4 * (1.742)
        g_mo4_1 = g_m2 + 4 * g_h2o - 4 * g_h2 + 5 * (1.507)
        g_mo4_3 = g_mo4_2 - 0.27
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
        # 2. dissolution to m+3
        dg[6] = g_m3 + g_ch0 - g_mn4 - 3 * u
        dg[7] = g_m3 + g_ch1 - g_mn4 - 2 * u - 0.5 * g_h2 + 1 * k * ph
        dg[8] = g_m3 + g_ch2 - g_mn4 - 1 * u - 1.0 * g_h2 + 2 * k * ph
        dg[9] = g_m3 + g_ch3 - g_mn4 - 0 * u - 1.5 * g_h2 + 3 * k * ph
        dg[10] = g_m3 + g_ch4 - g_mn4 + 1 * u - 2.0 * g_h2 + 4 * k * ph
        # 3. dissolution to g_mo4_2
        dg[11] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[12] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[13] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[14] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[15] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 5. dissolution to g_mo4_1
        dg[16] = (
            (g_mo4_1 + 4 * g_h2 - 8 * k * ph - 5 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[17] = (
            (g_mo4_1 + 4 * g_h2 - 8 * k * ph - 5 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[18] = (
            (g_mo4_1 + 4 * g_h2 - 8 * k * ph - 5 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[19] = (
            (g_mo4_1 + 4 * g_h2 - 8 * k * ph - 5 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[20] = (
            (g_mo4_1 + 4 * g_h2 - 8 * k * ph - 5 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 5. dissolution to g_mo4_3
        dg[21] = (
            ((g_mo4_3 + u) + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[22] = (
            ((g_mo4_3 + u) + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[23] = (
            ((g_mo4_3 + u) + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[24] = (
            ((g_mo4_3 + u) + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[25] = (
            ((g_mo4_3 + u) + 4 * g_h2 - 8 * k * ph - 4 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
    elif m == "Fe":
        g_m2 = g_m + 2 * (-0.447 + (1 / 2) * 0.0592 * np.log10(conc))
        g_m3 = g_m + 3 * (-0.037 + (1 / 3) * 0.0592 * np.log10(conc))
        g_hmo4_1 = g_m3 + 4 * g_h2o - (3 / 2) * g_h2 + 3 * (2.07)
        g_mo4_2 = g_m3 + 4 * g_h2o + 3 * (2.20) - 4 * g_h2
        g_moh2 = g_m3 + g_h2o - 0.5 * g_h2 + 2.43 * k
        g_hmo2_1 = g_m2 + 2 * g_h2o - (3 / 2) * g_h2 + 31.58 * k
        g_mo2_1 = g_hmo2_1 - (1 / 2) * g_h2 - 0.685
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
        # 2. dissolution to m+3
        dg[6] = g_m3 + g_ch0 - g_mn4 - 3 * u
        dg[7] = g_m3 + g_ch1 - g_mn4 - 2 * u - 0.5 * g_h2 + 1 * k * ph
        dg[8] = g_m3 + g_ch2 - g_mn4 - 1 * u - 1.0 * g_h2 + 2 * k * ph
        dg[9] = g_m3 + g_ch3 - g_mn4 - 0 * u - 1.5 * g_h2 + 3 * k * ph
        dg[10] = g_m3 + g_ch4 - g_mn4 + 1 * u - 2.0 * g_h2 + 4 * k * ph
        # 3. dissolution to hmo4_1
        dg[11] = (
            (g_hmo4_1 - 4 * g_h2o + (3 / 2) * g_h2 - 7 * k * ph - 3 * u)
            + g_ch0
            - g_mn4
            - 3 * u
        )
        dg[12] = (
            (g_hmo4_1 - 4 * g_h2o + (3 / 2) * g_h2 - 7 * k * ph - 3 * u)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[13] = (
            (g_hmo4_1 - 4 * g_h2o + (3 / 2) * g_h2 - 7 * k * ph - 3 * u)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[14] = (
            (g_hmo4_1 - 4 * g_h2o + (3 / 2) * g_h2 - 7 * k * ph - 3 * u)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[15] = (
            (g_hmo4_1 - 4 * g_h2o + (3 / 2) * g_h2 - 7 * k * ph - 3 * u)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 4. dissolution to mo4_2
        dg[16] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch0
            - g_mn4
            - 3 * u
        )
        dg[17] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[18] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[19] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[20] = (
            (g_mo4_2 + 4 * g_h2 - 8 * k * ph - 3 * u - 4 * g_h2o)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 5. dissolution to moh2
        dg[21] = (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph) + g_ch0 - g_mn4 - 3 * u
        dg[22] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch1
            - g_mn4
            - 2 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[23] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch2
            - g_mn4
            - 1 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[24] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch3
            - g_mn4
            - 0 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[25] = (
            (g_moh2 + 0.5 * g_h2 - g_h2o - k * ph)
            + g_ch4
            - g_mn4
            + 1 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 6. dissolution to hmo2_1
        dg[26] = (
            (g_hmo2_1 - 2 * g_h2o + (3 / 2) * g_h2 - 3 * k * ph)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[27] = (
            (g_hmo2_1 - 2 * g_h2o + (3 / 2) * g_h2 - 3 * k * ph)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[28] = (
            (g_hmo2_1 - 2 * g_h2o + (3 / 2) * g_h2 - 3 * k * ph)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[29] = (
            (g_hmo2_1 - 2 * g_h2o + (3 / 2) * g_h2 - 3 * k * ph)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[30] = (
            (g_hmo2_1 - 2 * g_h2o + (3 / 2) * g_h2 - 3 * k * ph)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
        # 7. dissolution to mo2_1
        dg[31] = (
            (
                (g_mo2_1 - (1 / 2) * g_h2 - 1 * k * ph)
                - 2 * g_h2o
                + (3 / 2) * g_h2
                - 3 * k * ph
            )
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[32] = (
            (
                (g_mo2_1 - (1 / 2) * g_h2 - 1 * k * ph)
                - 2 * g_h2o
                + (3 / 2) * g_h2
                - 3 * k * ph
            )
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[33] = (
            (
                (g_mo2_1 - (1 / 2) * g_h2 - 1 * k * ph)
                - 2 * g_h2o
                + (3 / 2) * g_h2
                - 3 * k * ph
            )
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[34] = (
            (
                (g_mo2_1 - (1 / 2) * g_h2 - 1 * k * ph)
                - 2 * g_h2o
                + (3 / 2) * g_h2
                - 3 * k * ph
            )
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[35] = (
            (
                (g_mo2_1 - (1 / 2) * g_h2 - 1 * k * ph)
                - 2 * g_h2o
                + (3 / 2) * g_h2
                - 3 * k * ph
            )
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
    elif m == "Co":
        g_m2 = g_m + 2 * (-0.28 + (1 / 2) * 0.0592 * np.log10(conc))
        g_m3 = g_m2 + 1.92
        g_hmo2_1 = g_m2 + 2 * g_h2o - (3 / 2) * g_h2 + 31.70 * k
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
        # 2. dissolution to m+3
        dg[6] = g_m3 + g_ch0 - g_mn4 - 3 * u
        dg[7] = g_m3 + g_ch1 - g_mn4 - 2 * u - 0.5 * g_h2 + 1 * k * ph
        dg[8] = g_m3 + g_ch2 - g_mn4 - 1 * u - 1.0 * g_h2 + 2 * k * ph
        dg[9] = g_m3 + g_ch3 - g_mn4 - 0 * u - 1.5 * g_h2 + 3 * k * ph
        dg[10] = g_m3 + g_ch4 - g_mn4 + 1 * u - 2.0 * g_h2 + 4 * k * ph
        # 3. dissolution to hmo2_1
        dg[11] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[12] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[13] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[14] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[15] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
    elif m == "Ni":
        g_m2 = g_m + 2 * (-0.26 + (1 / 2) * 0.0592 * np.log10(conc))
        g_m3 = g_m2 + 2.30
        g_hmo2_1 = g_m2 + 2 * g_h2o - (3 / 2) * g_h2 + 30.40 * k
        # 1. dissolution to m+2
        dg[1] = g_m2 + g_ch0 - g_mn4 - 2 * u - 0.0 * g_h2 + 0 * k * ph
        dg[2] = g_m2 + g_ch1 - g_mn4 - 1 * u - 0.5 * g_h2 + 1 * k * ph
        dg[3] = g_m2 + g_ch2 - g_mn4 - 0 * u - 1.0 * g_h2 + 2 * k * ph
        dg[4] = g_m2 + g_ch3 - g_mn4 + 1 * u - 1.5 * g_h2 + 3 * k * ph
        dg[5] = g_m2 + g_ch4 - g_mn4 + 2 * u - 2.0 * g_h2 + 4 * k * ph
        # 2. dissolution to m+3
        dg[6] = g_m3 + g_ch0 - g_mn4 - 3 * u
        dg[7] = g_m3 + g_ch1 - g_mn4 - 2 * u - 0.5 * g_h2 + 1 * k * ph
        dg[8] = g_m3 + g_ch2 - g_mn4 - 1 * u - 1.0 * g_h2 + 2 * k * ph
        dg[9] = g_m3 + g_ch3 - g_mn4 - 0 * u - 1.5 * g_h2 + 3 * k * ph
        dg[10] = g_m3 + g_ch4 - g_mn4 + 1 * u - 2.0 * g_h2 + 4 * k * ph
        # 3. dissolution to hmo2_1
        dg[11] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch0
            - g_mn4
            - 2 * u
            - 0.0 * g_h2
            + 0 * k * ph
        )
        dg[12] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch1
            - g_mn4
            - 1 * u
            - 0.5 * g_h2
            + 1 * k * ph
        )
        dg[13] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch2
            - g_mn4
            - 0 * u
            - 1.0 * g_h2
            + 2 * k * ph
        )
        dg[14] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch3
            - g_mn4
            + 1 * u
            - 1.5 * g_h2
            + 3 * k * ph
        )
        dg[15] = (
            (g_hmo2_1 + (3 / 2) * g_h2 - 3 * k * ph - 2 * g_h2o)
            + g_ch4
            - g_mn4
            + 2 * u
            - 2.0 * g_h2
            + 4 * k * ph
        )
    else:
        raise ValueError(f"Metal {m} not supported")
    dg[0] = 0  # bare surface
    dg[41] = g_mn4_ooh + (1.5 * g_h2) - 1 * k * ph - 1 * u - g_mn4 - (2 * g_h2o)

    ref = min(dg[1:41])  # minimum energy of ionic components
    b_rel = dg[0] - ref
    ooh_rel = dg[41] - ref
    return b_rel, ooh_rel


def get_master(data: dict) -> list:
    """Get the master xch database object and remaining energies from databases.

    Args:
        data (dict): Dictionary containing the data for the calculation.

    Returns:
        list: List of ASE database objects.
    """
    master_dir = Path(config["master_database_dir"])
    structure = data["name"]
    run_struc = data["run_structure"]
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    db_names = [
        "e_xch_solv_implicit_without_corrections",
        "e_xc_solv_implicit",
        "pristine_implicit",
        "e_adsorbate_without_corrections",
    ]
    mname = str(Path(run_struc).stem).replace(metal, "M")  # type: ignore

    xch_db = connect(master_dir / f"{db_names[0]}_master.db")

    xc_db = connect(master_dir / f"{db_names[1]}_master.db")
    e_xc = xc_db.get(name=mname).energy

    pristine_implicit = connect(master_dir / f"{db_names[2]}_master.db")
    e_mn4 = pristine_implicit.get(
        name=structure, metal=metal, carbon_structure=carbon_structure
    ).energy

    adsorption_db = connect(master_dir / f"{db_names[3]}_master.db")
    e_mn4_ooh = adsorption_db.get(name=structure, ads1="non", ads2="OOH").energy

    return [xch_db, e_xc, e_mn4, e_mn4_ooh]


def main(**data: dict) -> tuple[bool, dict | None]:
    """Calculate the relative stability of the OOH intermediate.

    Args:
        data: Dictionary containing the data for the calculation.

    Returns:
        Perqueue tuple containing boolean and dictionary.
    """
    carbon_structure = data["carbon_structure"]
    metal = str(data["metal"])
    base_dir = Path(str(data["base_dir"]))
    database_dir = Path(base_dir) / "runs" / "databases"
    master_dir = Path(config["master_database_dir"])
    merge(database_dir, master_dir)
    run_list = get_master(data)
    xch_db, e_xc, e_mn4, e_mn4_ooh = run_list
    run_struc = data["run_structure"]  # type: ignore
    structure = data["name"]
    dopant = data["dopant"]
    name = str(Path(run_struc).stem).replace(metal, "M")  # type: ignore
    h_master = {"0H": e_xc}
    for hs in range(1, 5):
        find_config_min = []
        # for conf in range(1, 5):
        for row in xch_db.select(carbon_structure=carbon_structure):
            try:
                xch_correction = retreive_sql_correction(
                    os.path.join(database_dir, "e_xch_vib_solv_implicit_corrections"),
                    {"name": f"{name}_{hs}H_config1"},
                )
            except NameError:
                xch_correction = retreive_sql_correction(
                    os.path.join(
                        master_dir, "e_xch_vib_solv_implicit_corrections_master"
                    ),
                    {"name": f"{name}_{hs}H_config1"},
                )

            if xch_correction is not None:
                corrected_energy = row.energy + float(xch_correction)
                find_config_min.append(corrected_energy)
            else:
                run_logger(
                    f"Could not find {name}_{hs}H_config1 in the thermal corrections database",
                    str(__file__),
                    "error",
                )
        if find_config_min:
            h_master[f"{hs}H"] = min(find_config_min)
        else:
            # initialize with a high value so if xh calculation(s) dont converge, it will be too large
            h_master[f"{hs}H"] = 9999
    if h_master["1H"] == 9999:
        run_logger(
            f"Could not find any 1H configurations for {name}", str(__file__), "error"
        )
        raise ValueError(f"Could not find any 1H configurations for {name}")
    thermal_ooh = float(
        retreive_sql_correction(
            os.path.join(database_dir, "e_ads_vib_corrections_master"), data
        )
    )
    e_mn4_ooh_corr = (
        e_mn4_ooh + thermal_ooh + 0.2
    )  # 0.2 eV correction for the Christensen correction
    b_rel, ooh_rel = acid_stability(ph, u, metal, e_mn4, h_master, e_mn4_ooh_corr)
    cutoff = 1.5
    print(f"b_rel: {b_rel}, ooh_rel: {ooh_rel}", flush=True)
    if b_rel < cutoff:
        return True, data
    else:
        run_logger(
            f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} -"
            + f"Operating stability b_rel; {b_rel} > {cutoff} eV. Ref: ooh_rel; {ooh_rel}",
            str(__file__),
            "error",
        )
        return False, None
