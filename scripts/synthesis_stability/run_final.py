import os
import re
from pathlib import Path

from ase.db import connect  # type: ignore
from utils import add_entry  # type: ignore


def main(**data: dict) -> tuple[bool, dict]:
    """Read all the individual terms and calculate the operating stability of the catalyst.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the base directory, run structure, carbon structure and metal.

    Returns:
        Perqueue return tuple.
    """
    # idx, *_ = data[INDEX_KW]
    # idx = str(idx)
    # data = data[idx]
    database = {}
    base_dir = data["base_dir"]
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    name = str(Path(str(data["run_structure"])).stem)
    data_base_folder = Path(str(base_dir)) / "runs" / "databases"
    os.system(
        f'cp {os.path.join(str(base_dir), "stock_databases", "e_m.db")}  {data_base_folder}'
    )
    e_m_on_c_db = connect(os.path.join(data_base_folder, "e_m_on_c.db"))
    e_c_db = connect(os.path.join(data_base_folder, "e_c.db"))
    e_mxc_db = connect(os.path.join(data_base_folder, "e_mxc.db"))
    e_m_db = connect(os.path.join(data_base_folder, "e_m.db"))
    e_xc_db = connect(os.path.join(data_base_folder, "e_xc.db"))
    e_xc = e_xc_db.get(name=name.replace(str(metal), "M")).energy
    e_c = e_c_db.get(carbon_structure=carbon_structure).energy
    e_mxc = e_mxc_db.get(
        metal=metal, carbon_structure=carbon_structure, name=name
    ).energy
    e_m_on_c = e_m_on_c_db.get(metal=metal, carbon_structure=carbon_structure).energy
    e_m = e_m_db.get(metal=metal).energy / int(
        re.findall(r"\d+\.?\d*", e_m_db.get(metal=metal).formula)[0]
    )  # ev/atom
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc

    database_data = {
        "g_a": g_a,
        "g_d": g_d,
        "carbon_structure": carbon_structure,
        "metal": metal,
        "dopant": data["dopant"],
        "run_structure": data["run_structure"],
        "base_dir": str(base_dir),
    }
    database[name] = database_data
    add_entry(os.path.join(data_base_folder, "operating_stability.json"), database)
    del data["pq_index"]

    # run_logger("DISCARD - Operating stability < ???? eV.", str(__file__), 'error'))
    # if 'z' in carbon_structure:
    #     add_entry(os.path.join(data_base_folder, "seperated_operating_stability.json"), database)
    #     return True, data
    # else:
    #     return False, None
    return True, data


if __name__ == "__main__":
    main()
