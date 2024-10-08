import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import (  # type: ignore
    check_ase_database,
    operating_stability_run_vasp,
    read_and_write_database,
    run_logger,
)
from vasp_input import vasp_input  # type: ignore


def e_xc(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_xc with solvation and dipole correction calculations.

    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        converged (bool): True if the calculation converged, raise error otherwise.
    """
    struc_path = Path(data["run_structure"])
    base_dir = Path(data["base_dir"])
    metal = data["metal"]
    del data["metal"]
    run_dir = Path(
        os.path.join(base_dir, "runs", "operating_stability", "e_xc_solv_implicit")
    )
    run_dir.mkdir(exist_ok=True, parents=True)
    structure = read(os.path.join(struc_path, "init.POSCAR"))
    del structure[[atom.symbol == metal for atom in structure]]
    remove_metal_dir = run_dir / Path(str(struc_path.stem).replace(metal, "M"))
    remove_metal_dir.mkdir(exist_ok=True, parents=True)
    structure.write(remove_metal_dir / "init.POSCAR")
    data["name"] = str(Path(data["run_structure"]).stem).replace(metal, "M")
    name = str(Path(data["run_structure"]).stem).replace(metal, "M")
    outcar = Path(remove_metal_dir) / "OUTCAR.RDip"
    if os.path.exists(remove_metal_dir / "OUTCAR.RDip"):
        print(
            f"e_xc_solv_implicit calculation already exists for {name} in {remove_metal_dir}",
            flush=True,
        )
        if not check_ase_database("e_xc_solv_implicit", data, master=False):
            print(f"Writing to local database for {name}", flush=True)
            read_and_write_database(outcar, "e_xc_solv_implicit", data)
        return True
    if check_ase_database("e_xc_solv_implicit", data, master=True):
        print("In master database already", flush=True)
        return True
    converged = operating_stability_run_vasp(
        remove_metal_dir, vasp_parameters, "e_xc_solv_implicit", data
    )
    if converged:
        read_and_write_database(outcar, "e_xc_solv_implicit", data)
        return True
    else:
        run_logger(
            "e_xc_solv_implicit calculation did not converge.", str(__file__), "error"
        )
        print("e_xc_solv_implicit calculation did not converge.", flush=True)
        raise ValueError("e_xc_implicit_solv calculation did not converge.")


def main(**data: dict) -> tuple[bool, None]:
    """Run the operating stability part of the workflow.

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
