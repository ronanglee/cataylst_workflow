import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import (  # type: ignore
    read_and_write_database,
    run_logger,
    synthesis_stability_run_vasp,
)
from vasp_input import vasp_input  # type: ignore


def e_c(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_c calculations.

    Args:
        data (dict): Dictionary containing the base directory and carbon structure.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        True if all the SCF calculations have converged, False otherwise.
    """
    carbon_structure = Path(data["carbon_structure"])
    base_dir = Path(data["base_dir"])
    e_c_dir = Path(
        os.path.join(
            base_dir, "runs", "synthesis_stability", "e_c", f"{carbon_structure}"
        )
    )
    e_c_dir.mkdir(parents=True, exist_ok=True)
    struct_dir = Path(
        os.path.join(base_dir, "template_structures", "carbon", f"{carbon_structure}")
    )
    structure = read(os.path.join(struct_dir, "init.POSCAR"))
    cwd = os.getcwd()
    structure.write(e_c_dir / "init.POSCAR")
    if os.path.exists(e_c_dir / "OUTCAR.opt"):
        outcar = Path(e_c_dir) / "OUTCAR.opt"
        data["name"] = data["carbon_structure"]
        return True
    converged = synthesis_stability_run_vasp(e_c_dir, vasp_parameters, "e_c")
    if converged:
        outcar = Path(e_c_dir) / "OUTCAR.opt"
        data["name"] = data["carbon_structure"]
        read_and_write_database(outcar, "e_c", data)
        os.chdir(cwd)
        return True
    else:
        print(f"e_c calculation did not converge for {carbon_structure}.")
        run_logger(f"e_c calculation did not converge for {carbon_structure}.", str(__file__), "error")
        os.chdir(cwd)
        raise ValueError(f"e_c calculation did not converge for {carbon_structure}.")


def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the run parameters.

    Returns:
        Perqueue return tuple.
    """
    vasp_parameters = vasp_input()
    converged = e_c(data, vasp_parameters)
    return converged, None


if __name__ == "__main__":
    main()
