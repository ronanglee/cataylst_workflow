import os
from pathlib import Path

from ase.db import connect  # type: ignore
from utils import read_config, run_logger  # type: ignore

config = read_config()


def try_local_then_master(data_base_folder: os.PathLike, data: dict) -> list:
    """Try to connect to the local database first, then the master database. If the local database
    does not exist, the master database will be used.

    Args:
        data_base_folder (os.PathLike): Path to the local database.
        data (dict): Dictionary containing the run parameters.

    Returns:
        list: List of energies.
    """
    metal = data["metal"]
    carbon_structure = data["carbon_structure"]
    name = str(Path(str(data["run_structure"])).stem)
    master_dir = config["master_database_dir"]
    local_dir = data_base_folder
    db_names = [
        "e_m_on_c",
        "e_c",
        "e_mxc",
        "e_m",
        "e_xc",
    ]
    try:
        e_m_on_c_db = connect(os.path.join(local_dir, f"{db_names[0]}.db"))
        e_m_on_c = e_m_on_c_db.get(
            metal=metal, carbon_structure=carbon_structure
        ).energy
    except AssertionError:
        e_m_on_c_db = connect(os.path.join(master_dir, f"{db_names[0]}_master.db"))
        e_m_on_c = e_m_on_c_db.get(
            metal=metal, carbon_structure=carbon_structure
        ).energy

    try:
        e_c_db = connect(os.path.join(local_dir, f"{db_names[1]}.db"))
        e_c = e_c_db.get(carbon_structure=carbon_structure).energy
    except AssertionError:
        e_c_db = connect(os.path.join(master_dir, f"{db_names[1]}_master.db"))
        e_c = e_c_db.get(carbon_structure=carbon_structure).energy

    try:
        e_mxc_db = connect(os.path.join(local_dir, f"{db_names[2]}.db"))
        e_mxc = e_mxc_db.get(
            metal=metal, carbon_structure=carbon_structure, name=name
        ).energy
    except AssertionError:
        e_mxc_db = connect(os.path.join(master_dir, f"{db_names[2]}_master.db"))
        e_mxc = e_mxc_db.get(
            metal=metal, carbon_structure=carbon_structure, name=name
        ).energy

    try:
        e_m_db = connect(os.path.join(local_dir, f"{db_names[3]}.db"))
        e_m = e_m_db.get(metal=metal).energy / len(
            e_m_db.get_atoms(metal=metal)
        )  # ev/atom
    except AssertionError:
        e_m_db = connect(os.path.join(master_dir, f"{db_names[3]}_master.db"))
        e_m = e_m_db.get(metal=metal).energy / len(
            e_m_db.get_atoms(metal=metal)
        )  # ev/atom

    try:
        e_xc_db = connect(os.path.join(local_dir, f"{db_names[4]}.db"))
        e_xc = e_xc_db.get(
            carbon_structure=carbon_structure, name=name.replace(str(metal), "M")
        ).energy
    except AssertionError:
        e_xc_db = connect(os.path.join(master_dir, f"{db_names[4]}_master.db"))
        e_xc = e_xc_db.get(
            carbon_structure=carbon_structure, name=name.replace(str(metal), "M")
        ).energy

    return [e_m_on_c, e_c, e_mxc, e_m, e_xc]


def main(**data: dict) -> tuple[bool, dict | None]:
    """Read all the individual terms and calculate the synthesis stability of the catalyst.
    g_a = -1 * (e_xc + e_m_on_c - e_c - e_mxc)
    hof = e_mxc - e_m - e_xc
    where x is the dopant, c is carbon and m is the metal.

    g_a - approximated intermediate structure energy
    g_d - approximated final structure energy
    hof - heat of formation

    Args:
        data (dict): Dictionary containing the run parameters.

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
    vals = try_local_then_master(data_base_folder, data)
    e_m_on_c, e_c, e_mxc, e_m, e_xc = vals
    g_a = -1 * (e_xc + e_m_on_c - e_c - e_mxc)
    hof = e_mxc - e_m - e_xc
    print(f"g_a: {g_a}, hof: {hof}", flush=True)
    if hof > 0:
        run_logger(
            f"DISCARD - {structure}/{carbon_structure}/{dopant}/{metal} - Heat of formation; {hof} > 0 eV",
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
