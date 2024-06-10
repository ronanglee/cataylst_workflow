import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import read_and_write_datbase  # type: ignore
from vasp_input import synthesis_stability_run_vasp, vasp_input  # type: ignore


def e_mxc(data: dict, vasp_parameters: dict) -> None:
    """Generate and run the input files for the e_mxc calculations.

    Args:
        data (dict): Dictionary containing the run structure and base directory.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    struc_path = Path(data["run_structure"])
    base_dir = Path(data["base_dir"])
    structure_name = Path(struc_path).stem
    e_mxc_dir = Path(
        os.path.join(
            base_dir, "runs", "synthesis_stability", "e_mxc", f"{structure_name}"
        )
    )
    e_mxc_dir.mkdir(parents=True, exist_ok=True)
    structure = read(os.path.join(struc_path, "init.POSCAR"))
    structure.write(e_mxc_dir / "init.POSCAR")
    cwd = os.getcwd()
    synthesis_stability_run_vasp(e_mxc_dir, vasp_parameters)
    read_and_write_datbase(e_mxc_dir, base_dir, "e_mxc", data)

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
    e_mxc(data, vasp_parameters)

    return True, None


if __name__ == "__main__":
    main()
