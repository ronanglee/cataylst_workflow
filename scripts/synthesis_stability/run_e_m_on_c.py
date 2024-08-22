import os
from pathlib import Path

from ase.io import read  # type: ignore
from utils import (  # type: ignore
    gather_structs,
    read_and_write_database,
    run_logger,
    synthesis_stability_run_vasp,
    check_database
)
from vasp_input import vasp_input  # type: ignore


def e_m_on_c(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_m_at_c calculations.

    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        True if all the SCF calculations have converged, False otherwise.
    """
    skimmed_data = data.copy()
    del skimmed_data["all_run_structures"]
    del skimmed_data["dopants"]
    carbon_structure = Path(skimmed_data["carbon_structure"])
    base_dir = Path(skimmed_data["base_dir"])
    metal = skimmed_data["metal"]
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
    skimmed_data["name"] = f"{metal}{carbon_structure}_0N_0H"
    if os.path.exists(e_m_on_c_dir / "OUTCAR.opt"):
        outcar = Path(e_m_on_c_dir) / "OUTCAR.opt"
        run_logger(f"e_m_on_c calculation already exists for {metal} on {carbon_structure}.", str(__file__), "info")
        print(f"e_m_on_c calculation already exists for {metal} on {carbon_structure}.")
        return True
    if check_database("e_m_on_c", skimmed_data):
        print('In master database already')        
        return True
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
    structure.write(e_m_on_c_dir / "init.POSCAR")
    converged = synthesis_stability_run_vasp(e_m_on_c_dir, vasp_parameters, "e_m_on_c")
    if converged:
        outcar = Path(e_m_on_c_dir) / "OUTCAR.opt"
        skimmed_data["name"] = f"{metal}_{carbon_structure}_0N_0H"
        read_and_write_database(outcar, "e_m_on_c", skimmed_data)
        os.chdir(cwd)
        return True
    else:
        run_logger(f"e_m_on_c calculation did not converge for {metal} on {carbon_structure}.", str(__file__), "error")
        print(f"e_m_on_c calculation did not converge for {metal} on {carbon_structure}.")
        os.chdir(cwd)
        raise ValueError(f"e_m_on_c calculation did not converge for {metal} on {carbon_structure}.")


def main(**data: dict) -> tuple[bool, dict]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing the following keys:
            base_dir (str): Path to the base workflow directory.
            run_structure (str): Path to the generated input files.
            carbon_structure str): Identity of carbon structure (bulk, armchair or zigzag).
            metal (str): Metal in structure.

    Returns:
        Perqueue return tuple.
    """
    vasp_parameters = vasp_input()
    metals = data["metals"]
    del data["metals"]
    converged = e_m_on_c(data, vasp_parameters)
    workflow_data = gather_structs(data, metals)
    return converged, workflow_data


if __name__ == "__main__":
    main()
