import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import read_and_write_datbase  # type: ignore
from vasp_input import synthesis_stability_run_vasp, vasp_input  # type: ignore


def e_xc(
    struc_path: os.PathLike, base_dir: os.PathLike, metal: str, vasp_parameters: dict
) -> None:
    """Generate and run the input files for the e_xc calculations.

    Args:
        strucs_path (Path): Path to the generated input file.
        base_dir (Path): Path to the workflow script directory.
        metal (str): Metal in structure.
        vas_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    run_dir = Path(os.path.join(base_dir, "runs", "synthesis_stability", "e_xc"))
    run_dir.mkdir(exist_ok=True, parents=True)
    structure = read(os.path.join(struc_path, "init.POSCAR"))
    del structure[[atom.symbol == metal for atom in structure]]
    remove_metal_dir = run_dir / struc_path / "remove_metal"
    remove_metal_dir.mkdir(exist_ok=True, parents=True)
    structure.write(remove_metal_dir / "init.POSCAR")
    synthesis_stability_run_vasp(remove_metal_dir, vasp_parameters)
    read_and_write_datbase(remove_metal_dir, base_dir)
    os.chdir(base_dir)


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
            dopant (List): List of dopants. ???? Is this needed?

    Returns:
        Perqueue tuple containing a boolean and None.
    """
    vasp_parameters = vasp_input()
    e_xc(Path(data["run_structure"]), Path(data["base_dir"]), data["metal"], vasp_parameters)  # type: ignore

    return True, None


if __name__ == "__main__":
    main()
