import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import read_and_write_datbase  # type: ignore
from vasp_input import synthesis_stability_run_vasp, vasp_input  # type: ignore


def e_c(data: dict, vasp_parameters: dict) -> None:
    """Generate and run the input files for the e_c calculations.

    Args:
        data (dict): Dictionary containing the base directory and carbon structure.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
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
    synthesis_stability_run_vasp(e_c_dir, vasp_parameters)
    read_and_write_datbase(e_c_dir, base_dir, "e_c")
    os.chdir(cwd)


def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the following keys:
            base_dir (str): Path to the base workflow directory.
            run_structure (str): Path to the generated input files.
            carbon_structure (str): Identity of carbon structure (bulk, armchair or zigzag).
            metals (str): Metal in structure.

    Returns:
        Perqueue tuple containing a boolean and None.
    """
    vasp_parameters = vasp_input()
    e_c(data, vasp_parameters)

    return True, None


if __name__ == "__main__":
    main()
