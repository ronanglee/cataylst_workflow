import os
from pathlib import Path

from ase.io import read  # type: ignore
from perqueue.constants import INDEX_KW  # type: ignore
from utils import (  # type: ignore
    read_and_write_database,
    run_logger,
    operating_stability_run_vasp,
)
from vasp_input import vasp_input  # type: ignore


def e_xc(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_xc with solvation and dipole correction calculations.

    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        bool: True if the calculation converged, False otherwise.
    """
    struc_path = Path(data["run_structure"])
    base_dir = Path(data["base_dir"])
    metal = data["metal"]
    dopant = data["dopant"]
    del data["metal"]
    run_dir = Path(os.path.join(base_dir, "runs", "operating_stability", "e_xc_solv_implicit"))
    run_dir.mkdir(exist_ok=True, parents=True)
    structure = read(os.path.join(struc_path, "init.POSCAR"))
    if dopant != "":
        if dopant != 'O':
            if dopant != 'SB':
                for atom in structure:
                    if atom.symbol == "B":
                        atom.symbol = dopant
    del structure[[atom.symbol == metal for atom in structure]]
    remove_metal_dir = run_dir / Path(str(struc_path.stem).replace(metal, "M"))
    remove_metal_dir.mkdir(exist_ok=True, parents=True)
    structure.write(remove_metal_dir / "init.POSCAR")
    if os.path.exists(remove_metal_dir / "OUTCAR.RDip"):
        outcar = Path(remove_metal_dir) / "OUTCAR.RDip"
        data["name"] = str(Path(data["run_structure"]).stem).replace(metal, "M")
        return True
    converged = operating_stability_run_vasp(remove_metal_dir, vasp_parameters, "e_xc_solv_implicit",)
    if converged:
        outcar = Path(remove_metal_dir) / "OUTCAR.RDip"
        data["name"] = str(Path(data["run_structure"]).stem).replace(metal, "M")
        read_and_write_database(outcar, "e_xc_solv_implicit", data)
        return True
    else:
        run_logger("e_xc_solv_implicit calculation did not converge.", str(__file__), "error")
        print("e_xc_solv_implicit calculation did not converge.")
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
