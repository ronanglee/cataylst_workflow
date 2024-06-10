import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import read_and_write_datbase  # type: ignore
from vasp_input import synthesis_stability_run_vasp, vasp_input  # type: ignore


def e_c(
    carbon_structure: os.PathLike, base_dir: os.PathLike, vasp_parameters: dict
) -> None:
    """Generate and run the input files for the e_c calculations.

    Args:
        carbon_structure (Path): Identity of carbon structure (bulk, armchair or zigzag).
        base_dir (Path): Path to the workflow script directory.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
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
    read_and_write_datbase(e_c_dir, base_dir)
    os.chdir(cwd)


def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mnc
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
    e_c(Path(data["carbon_structure"]), Path(data["base_dir"]), vasp_parameters)  # type: ignore

    return True, None


if __name__ == "__main__":
    main()
