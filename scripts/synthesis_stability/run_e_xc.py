from vasp_input import synthesis_stability_run_vasp
from vasp_input import vasp_input
import os
import sys
import posixpath
from ase.calculators.vasp import Vasp

from pathlib import Path
from ase.io import read

sys.path.append(str(Path(__file__).parent.parent))


def e_xc(
    struc_path: str, base_dir: os.PathLike, metal: str, vasp_parameters: dict
) -> None:
    """Generate and run the input files for the e_xc calculations.

    Args:
        strucs_path (str): Path to the generated input file.
        base_dir (Path): Path to the workflow script directory.
        metal (Str): Metal in structure.
        vas_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    run_dir = base_dir / "runs" / "synthesis_stability" / "e_xc"
    run_dir.mkdir(exist_ok=True, parents=True)
    struc_dir = Path(struc_path)
    structure = read(struc_dir / "init.POSCAR")
    del structure[[atom.symbol == metal for atom in structure]]
    remove_metal_dir = run_dir / posixpath.basename(struc_dir) / "remove_metal"
    remove_metal_dir.mkdir(exist_ok=True, parents=True)
    structure.write(remove_metal_dir / "init.POSCAR")
    synthesis_stability_run_vasp(remove_metal_dir, vasp_parameters)
    os.chdir(base_dir)


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
    run_dir = base_dir / "runs" / "synthesis_stability" / "e_m_on_c"
    e_m_on_c_dir = run_dir / f"{carbon_structure}_0N_0H" / f"{metal}"
    if e_m_on_c_dir.exists():  # Skip as carbon-metal structure already exists
        return
    else:
        e_m_on_c_dir.mkdir(exist_ok=True, parents=True)
    struct_dir = (
        base_dir / "template_structures" / "m_on_c" / f"{carbon_structure}_0N_0H"
    )
    structure = read(struct_dir / "init.POSCAR")
    for atom in structure:
        if atom.symbol == "Pt":  # Pt is the template metal
            atom.symbol = metal
    cwd = os.getcwd()
    print(f"creating {e_m_on_c_dir}")
    structure.write(e_m_on_c_dir / "init.POSCAR")
    run_vasp(e_m_on_c_dir, vasp_parameters)
    os.chdir(cwd)


def e_c(carbon_structure: str, base_dir: os.PathLike, vasp_parameters: dict) -> None:
    """Generate and run the input files for the e_c calculations.

    Args:
        carbon_structure (Str): Identity of carbon structure (bulk, armchair or zigzag).
        base_dir (Path): Path to the workflow script directory.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    run_dir = base_dir / "runs" / "synthesis_stability" / "e_c"
    e_c_dir = run_dir / f"{carbon_structure}"
    if e_c_dir.exists():  # Skip as carbon structure already exists
        return
    else:
        e_c_dir.mkdir(exist_ok=True, parents=True)
    struct_dir = base_dir / "template_structures" / "carbon" / f"{carbon_structure}"
    structure = read(struct_dir / "init.POSCAR")
    cwd = os.getcwd()
    print(f"creating {e_c_dir}")
    structure.write(e_c_dir / "init.POSCAR")
    run_vasp(e_c_dir, vasp_parameters)
    os.chdir(cwd)


def e_mxc(base_dir: os.PathLike, vasp_parameters: dict) -> None:
    """Generate and run the input files for the e_mxc calculations.

    Args:
        base_dir (Path): Path to the workflow script directory.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    run_dir = base_dir / "runs" / "synthesis_stability" / "e_mxc"
    e_mxc_dir = run_dir / f"{base_dir.name}"
    structure = read(base_dir / "init.POSCAR")
    structure.write(e_mxc_dir / "init.POSCAR")
    cwd = os.getcwd()
    print(f"creating {e_mxc_dir}")
    run_vasp(e_mxc_dir, vasp_parameters)
    os.chdir(cwd)


def main(**data) -> None:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mnc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the following keys:
            base_dir (Str): Path to the base workflow directory.
            run_structure (Str): Path to the generated input files.
            carbon_structure (Str): Identity of carbon structure (bulk, armchair or zigzag).
            metals (Str): Metal in structure.
            dopant (List): List of dopants. ???? Is this needed?

    Returns:
        None
    """
    vasp_parameters = vasp_input()
    e_xc(data["run_structure"], Path(data["base_dir"]), data["metal"], vasp_parameters)
    e_m_on_c(
        data["carbon_structure"], Path(data["base_dir"]), data["metal"], vasp_parameters
    )
    e_c(data["carbon_structure"], Path(data["base_dir"]), vasp_parameters)

    e_mxc(Path(data["run_structure"]), vasp_parameters)


if __name__ == "__main__":
    main(**data)
