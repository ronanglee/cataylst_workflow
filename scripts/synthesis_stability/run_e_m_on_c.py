import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import read_and_write_datbase  # type: ignore
from vasp_input import synthesis_stability_run_vasp, vasp_input  # type: ignore


def e_m_on_c(
    carbon_structure: str, base_dir: os.PathLike, metal: str, vasp_parameters: dict
) -> None:
    """Generate and run the input files for the e_m_at_c calculations.

    Args:
        carbon_structure (Str): Identity of carbon structure (bulk, armchair or zigzag).
        base_dir (Path): Path to the workflow script directory.
        metal (Str): Metal in structure.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    e_m_on_c_dir = Path(
        os.path.join(
            base_dir,
            "runs",
            "synthesis_stability",
            "e_m_on_c",
            f"{carbon_structure}_0N_0H",
            f"{metal}",
        )
    )
    e_m_on_c_dir.mkdir(parents=True, exist_ok=True)
    struct_dir = Path(
        os.path.join(
            base_dir, "template_structures", "m_on_c", f"{carbon_structure}_0N_0H"
        )
    )
    structure = read(struct_dir / "init.POSCAR")
    for atom in structure:
        if atom.symbol == "Pt":  # Pt is the template metal
            atom.symbol = metal
    cwd = os.getcwd()
    print(f"creating {e_m_on_c_dir}")
    structure.write(e_m_on_c_dir / "init.POSCAR")
    synthesis_stability_run_vasp(e_m_on_c_dir, vasp_parameters)
    read_and_write_datbase(e_m_on_c_dir, base_dir)
    os.chdir(cwd)


def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mnc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the following keys:
            base_dir (Path): Path to the base workflow directory.
            run_structure (Path): Path to the generated input files.
            carbon_structure (Str): Identity of carbon structure (bulk, armchair or zigzag).
            metals (Str): Metal in structure.

    Returns:
        Perqueue tuple containing a boolean and None.
    """
    vasp_parameters = vasp_input()
    e_m_on_c(  # type: ignore
        data["carbon_structure"], data["base_dir"], data["metal"], vasp_parameters  # type: ignore
    )  # type: ignore
    return True, None


if __name__ == "__main__":
    main()
