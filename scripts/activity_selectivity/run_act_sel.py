import os
from pathlib import Path  # type: ignore

from ase.db import connect  # type: ignore
from utils import read_config, retreive_sql_correction, run_logger  # type: ignore

g_h2 = -7.18 + 0.09  # eV/molecule (Free energy + BEEF correction)
g_h2o = -12.81 + (-0.03)  # eV/molecule (Free energy + BEEF correction)
dg_h2o = -4.92

database_dir = Path(__file__).parent.parent.parent / "runs" / "databases"

config = read_config()


def gibbs_ooh(g_ooh: float, g_bare: float, thermal_ooh: float) -> float:
    """Calculate the Gibbs free energy of OOH with 0.02 being
    the Christensen correction.

    Args:
        g_ooh (float): Free energy of OOH
        g_bare (float): Free energy of bare surface
        thermal_ooh (float): Thermal correction of OOH

    Returns:
        float: Gibbs free energy of OOH
    """
    return g_ooh + (1.5 * g_h2) - g_bare - (2 * g_h2o) + thermal_ooh + 0.2


def main(**data: dict) -> tuple[bool, dict | None]:
    """Calculate selectivity and activity

    Args:
        data (dict): Dictionary containing the run parameters.

    Returns:
        Perqueue return tuple
    """
    master_dir = Path(config["master_database_dir"])
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    dopant = data["dopant"]
    adsorption_db = connect(master_dir / "e_adsorbate_without_corrections_master.db")
    pristine_db = connect(master_dir / "pristine_implicit_master.db")
    thermal_ooh = float(
        retreive_sql_correction(
            os.path.join(master_dir, "e_ads_vib_corrections_master"), data
        )
    )
    structure = str(Path(data["run_structure"]).stem)  # type: ignore
    g_ooh = adsorption_db.get(name=structure, ads1="non", ads2="OOH").energy
    g_bare = pristine_db.get(name=structure).energy
    delta_g1 = gibbs_ooh(g_ooh, g_bare, thermal_ooh) + dg_h2o
    # 4.92 dg_h2o, d_g (delta_g) 1.4 h2o2 to water
    delta_g2 = (4.92 - 1.4) - gibbs_ooh(g_ooh, g_bare, thermal_ooh)
    ooh_descriptor = gibbs_ooh(g_ooh, g_bare, thermal_ooh)

    # Limiting potential (-1 to make it a volcano)
    ul_2e = -1 * max(delta_g1, delta_g2)
    del data["pq_index"]
    ul_2e_cutoff = 0.2
    ooh_descriptor_cutoffs = [3.5, 5]
    if (
        ooh_descriptor_cutoffs[0] < ooh_descriptor < ooh_descriptor_cutoffs[1]
        and ul_2e > ul_2e_cutoff
    ):
        return True, None
    else:
        if ul_2e < ul_2e_cutoff:
            if (
                ooh_descriptor < ooh_descriptor_cutoffs[0]
                or ooh_descriptor > ooh_descriptor_cutoffs[1]
            ):
                run_logger(
                    f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal} - activity ul_2e"
                    + f"{ul_2e} < {ul_2e_cutoff} eV and selectivity {ooh_descriptor} outside"
                    + f"[{ooh_descriptor_cutoffs[0]}, {ooh_descriptor_cutoffs[1]}] eV.",
                    str(__file__),
                    "error",
                )
            else:
                run_logger(
                    f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal}"
                    + f"- activity ul_2e {ul_2e} < {ul_2e_cutoff} eV",
                    str(__file__),
                    "error",
                )
        else:
            run_logger(
                f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal} -"
                + f"selectivity {ooh_descriptor} < {ooh_descriptor_cutoffs[0]}"
                + f"or > {ooh_descriptor_cutoffs[1]} eV.",
                str(__file__),
                "error",
            )
        return False, None
