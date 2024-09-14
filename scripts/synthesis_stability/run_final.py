import os
import re
from pathlib import Path

from ase.db import connect  # type: ignore
from utils import run_logger  # type: ignore


def main(**data: dict) -> tuple[bool, dict | None]:
    """Read all the individual terms and calculate the operating stability of the catalyst.
    g_a = -1 * (e_xc + e_m_on_c - e_c - e_mxc)
    hof = e_mxc - e_m - e_xc
    where x is the dopant, c is carbon and m is the metal.

    g_a - approximated intermediate structure energy
    g_d - approximated final structure energy
    hof - heat of formation

    Args:
        data (dict): Dictionary containing the base directory, run structure, carbon structure and metal.

    Returns:
        Perqueue return tuple.
    """
    base_dir = data["base_dir"]
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    structure = str(Path(str(data["run_structure"])).stem)
    dopant = data["dopant"]
    data_base_folder = Path(str(base_dir)) / "runs" / "databases"
    os.system(
        f'cp {os.path.join(str(base_dir), "stock_databases", "e_m.db")}  {data_base_folder}'
    )
    e_m_on_c_db = connect(os.path.join(data_base_folder, "e_m_on_c.db"))
    e_m_on_c = e_m_on_c_db.get(metal=metal, carbon_structure=carbon_structure).energy
    e_c_db = connect(os.path.join(data_base_folder, "e_c.db"))
    e_c = e_c_db.get(carbon_structure=carbon_structure).energy
    e_mxc_db = connect(os.path.join(data_base_folder, "e_mxc.db"))
    e_mxc = e_mxc_db.get(
        metal=metal, carbon_structure=carbon_structure, name=structure
    ).energy
    e_m_db = connect(os.path.join(data_base_folder, "e_m.db"))
    e_m = e_m_db.get(metal=metal).energy / int(
        re.findall(r"\d+\.?\d*", e_m_db.get(metal=metal).formula)[0]
    )  # ev/atom
    e_xc_db = connect(os.path.join(data_base_folder, "e_xc.db"))
    e_xc = e_xc_db.get(name=structure.replace(str(metal), "M")).energy
    g_a = -1 * (e_xc + e_m_on_c - e_c - e_mxc)
    hof = e_mxc - e_m - e_xc
    if hof > 0:
        run_logger(
            f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} - Heat of formation; {hof} < 0 eV",
            str(__file__),
            "error",
        )
        return False, None
    if g_a < 0:
        return True, data
    else:
        run_logger(
            f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} - Synthesis stability g_a; {g_a} > 0 eV.",
            str(__file__),
            "error",
        )
        return False, None


if __name__ == "__main__":
    main()
