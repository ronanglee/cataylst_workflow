import json
import os
from pathlib import Path  # type: ignore

from ase.db import connect  # type: ignore
from utils import add_entry, run_logger  # type: ignore

g_h2 = -7.18 + 0.09  # eV/molecule (Free energy + BEEF correction)
g_h2o = -12.81 + (-0.03)  # eV/molecule (Free energy + BEEF correction)
dg_h2o = -4.92

database_dir = Path(__file__).parent.parent.parent / "runs" / "databases"


def main(**data: dict) -> tuple[bool, dict | None]:
    """Calculate selectivity and activity

    Args:
        data (dict): Dictionary containing the run parameters.

    Returns:
        Perqueue return tuple
    """
    database = {}
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    dopant = data["dopant"]
    adsorption_db = connect(database_dir / "adsorption.db")
    pristine_db = connect(database_dir / "pristine_implicit.db")
    thermal_corrections = json.load(
        open(os.path.join(database_dir, "ads_vib_corrections.json"))
    )
    structure = str(Path(str(data["run_structure"])).stem)
    g_ooh = adsorption_db.get(name=structure, ads1="non", ads2="OOH").energy
    g_bare = pristine_db.get(name=structure).energy
    thermal_ooh = thermal_corrections[structure]["correction"]

    def gibbs_ooh(g_ooh: float, g_bare: float) -> float:
        """Calculate the Gibbs free energy of OOH with 0.02 being
        the Christensen correction.

        Args:
            g_ooh (float): Free energy of OOH
            g_bare (float): Free energy of bare surface

        Returns:
            float: Gibbs free energy of OOH
        """
        return g_ooh + (1.5 * g_h2) - g_bare - (2 * g_h2o) + thermal_ooh + 0.2

    delta_g1 = gibbs_ooh(g_ooh, g_bare) + dg_h2o
    # 4.92 dg_h2o, d_g (delta_g) 1.4 h2o2 to water
    delta_g2 = (4.92 - 1.4) - gibbs_ooh(g_ooh, g_bare)
    ooh_descriptor = gibbs_ooh(g_ooh, g_bare)

    # Limiting potential (-1 to make it a volcano)
    ul_2e = -1 * max(delta_g1, delta_g2)

    database_data = {
        "ooh_descriptor": ooh_descriptor,
        "ul_2e": ul_2e,
        "carbon_structure": carbon_structure,
        "metal": metal,
        "dopant": data["dopant"],
        "run_structure": data["run_structure"],
        "base_dir": str(data["base_dir"]),
    }
    database[structure] = database_data
    add_entry(os.path.join(database_dir, "act_sel.json"), database)
    del data["pq_index"]
    ul_2e_cutoff = 0.2
    ooh_descriptor_cutoffs = [3.5, 5]
    if (
        ooh_descriptor_cutoffs[0] < ooh_descriptor < ooh_descriptor_cutoffs[1]
        and ul_2e > ul_2e_cutoff
    ):
        add_entry(os.path.join(database_dir, "seperated_act_sel.json"), database)
        return True, None
    else:
        if ul_2e < ul_2e_cutoff:
            if (
                ooh_descriptor < ooh_descriptor_cutoffs[0]
                or ooh_descriptor > ooh_descriptor_cutoffs[1]
            ):
                string1 = f"selectivity {ooh_descriptor} outside  [{ooh_descriptor_cutoffs[0]}, {ooh_descriptor_cutoffs[1]}] eV."  # noqa
                run_logger(
                    f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal} - activity ul_2e {ul_2e} < {ul_2e_cutoff} eV and {string1}",  # noqa
                    str(__file__),
                    "error",
                )
            else:
                run_logger(
                    f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal} - activity ul_2e {ul_2e} < {ul_2e_cutoff} eV",  # noqa
                    str(__file__),
                    "error",
                )
        else:
            string1 = f"selectivity {ooh_descriptor} < {ooh_descriptor_cutoffs[0]} or > {ooh_descriptor_cutoffs[1]} eV."  # noqa
            run_logger(
                f"DISCARD - {structure}/{carbon_structure}/{dopant}{metal} - {string1}",
                str(__file__),
                "error",
            )
        return False, None
