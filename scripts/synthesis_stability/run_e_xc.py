import os
from pathlib import Path

from ase.io import read  # type: ignore
from perqueue.constants import INDEX_KW  # type: ignore
from utils import (  # type: ignore
    read_and_write_database,
    run_logger,
    synthesis_stability_run_vasp,
    check_database
)
from vasp_input import vasp_input  # type: ignore


def e_xc(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_xc calculations.

    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        bool: True if the calculation converged, False otherwise.
    """
    struc_path = Path(data["run_structure"])
    base_dir = Path(data["base_dir"])
    metal = data["metal"]
    del data["metal"]
    run_dir = Path(os.path.join(base_dir, "runs", "synthesis_stability", "e_xc"))
    run_dir.mkdir(exist_ok=True, parents=True)
    structure = read(os.path.join(struc_path, "init.POSCAR"))
    del structure[[atom.symbol == metal for atom in structure]]
    remove_metal_dir = run_dir / Path(str(struc_path.stem).replace(metal, "M"))
    remove_metal_dir.mkdir(exist_ok=True, parents=True)
    structure.write(remove_metal_dir / "init.POSCAR")
    data["name"] = str(Path(data["run_structure"]).stem).replace(metal, "M")
    if os.path.exists(remove_metal_dir / "OUTCAR.opt"):
        return True
    if check_database("e_xc", data):
        print('In master database already')
        return True
    converged = synthesis_stability_run_vasp(remove_metal_dir, vasp_parameters, "e_xc")
    if converged:
        outcar = Path(remove_metal_dir) / "OUTCAR.opt"
        read_and_write_database(outcar, "e_xc", data)
        return True
    else:
        run_logger("e_xc calculation did not converge.", str(__file__), "error")
        print("e_xc calculation did not converge.")
        raise ValueError("e_xc calculation did not converge.")

def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing run parameters.
    Returns:
        Perqueue return tuple.
    """
    cwd = os.getcwd()
    vasp_parameters = vasp_input()
    converged = e_xc(data, vasp_parameters)  # type: ignore
    os.chdir(cwd)
    return converged, None


if __name__ == "__main__":
    main()
