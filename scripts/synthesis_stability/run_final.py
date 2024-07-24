import os
import re
from pathlib import Path

from ase.db import connect  # type: ignore
from utils import add_entry, run_logger  # type: ignore


def main(**data: dict) -> tuple[bool, dict | None]:
    """Read all the individual terms and calculate the operating stability of the catalyst.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
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
    database = {}
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
    e_c_db = connect(os.path.join(data_base_folder, "e_c.db"))
    e_mxc_db = connect(os.path.join(data_base_folder, "e_mxc.db"))
    e_m_db = connect(os.path.join(data_base_folder, "e_m.db"))
    e_xc_db = connect(os.path.join(data_base_folder, "e_xc.db"))
    e_xc = e_xc_db.get(name=structure.replace(str(metal), "M")).energy
    e_c = e_c_db.get(carbon_structure=carbon_structure).energy
    e_mxc = e_mxc_db.get(
        metal=metal, carbon_structure=carbon_structure, name=structure
    ).energy
    e_m_on_c = e_m_on_c_db.get(metal=metal, carbon_structure=carbon_structure).energy
    e_m = e_m_db.get(metal=metal).energy / int(
        re.findall(r"\d+\.?\d*", e_m_db.get(metal=metal).formula)[0]
    )  # ev/atom
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    hof = e_mxc - e_m - e_xc
    database_data = {
        "g_a": g_a,
        "g_d": g_d,
        "hof": hof,
        "carbon_structure": carbon_structure,
        "metal": metal,
        "dopant": data["dopant"],
        "run_structure": data["run_structure"],
        "base_dir": str(base_dir),
    }
    database[structure] = database_data
    add_entry(os.path.join(data_base_folder, "synthesis_stability.json"), database)
    del data["pq_index"]

    if hof > 0:
        run_logger(f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} - Heat of formation; {hof} < 0 eV", str(__file__), 'error')
        return False, None
    if g_a > 0:
        add_entry(os.path.join(data_base_folder, "seperated_synthesis_stability.json"), database)
        return True, data
    else:
        run_logger(f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} - Synthesis stability g_a; {g_a} > 0 eV.", str(__file__), 'error')
        return False, None

if __name__ == "__main__":
    main()
