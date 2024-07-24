import os
from pathlib import Path

from ase.db import connect  # type: ignore
from ase.io import read  # type: ignore
from typing import Optional
from perqueue.constants import INDEX_KW  # type: ignore
from utils import (  # type: ignore
    read_and_write_database,
    run_logger,
    synthesis_stability_run_vasp,
)
from vasp_input import vasp_input  # type: ignore


def e_mxc(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_mxc calculations.

    Args:
        data (dict): Dictionary containing the run structure and base directory.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
       True if all the SCF calculations have converged, False otherwise.
    """
    struc_path = Path(data["run_structure"])
    base_dir = Path(data["base_dir"])
    dopant = data["dopant"]
    structure_name = Path(struc_path).stem
    dopant_template = Path(struc_path).stem.replace("B", dopant)
    dopant_template_path = struc_path.parent / dopant_template
    structure = read(os.path.join(dopant_template_path, "init.POSCAR"))
    e_mxc_dir = Path(
        os.path.join(
            base_dir, "runs", "synthesis_stability", "e_mxc", f"{structure_name}"
        )
    )
    if os.path.exists(e_mxc_dir):
        outcar = Path(e_mxc_dir) / "OUTCAR.opt"
        data["name"] = str(Path(data["run_structure"]).stem)
        return True
    e_mxc_dir.mkdir(parents=True, exist_ok=True)
    if dopant != "":
        if dopant != 'O':
            if dopant != 'SB':
                for atoms in structure:
                    if atoms.symbol == "B":
                        atoms.symbol = dopant
    structure.write(e_mxc_dir / "init.POSCAR")
    cwd = os.getcwd()
    converged = synthesis_stability_run_vasp(e_mxc_dir, vasp_parameters, "e_mxc")
    if converged:
        outcar = Path(e_mxc_dir) / "OUTCAR.opt"
        data["name"] = str(Path(data["run_structure"]).stem)
        read_and_write_database(outcar, "e_mxc", data)
        os.chdir(cwd)
        return True
    else:
        run_logger(f"e_mxc calculation did not converge for {structure}.", str(__file__), "error")
        print(f"e_mxc calculation did not converge for {structure}.")
        os.chdir(cwd)
        raise ValueError(f"e_mxc calculation did not converge for {structure}.")


def main(**data: dict) -> tuple[bool, Optional[dict]]:
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
    idx, *_ = data[INDEX_KW]
    idx = str(idx)
    data_base_folder = Path(data[idx]["base_dir"]) / "runs" / "databases"
    e_xc_db = connect(os.path.join(data_base_folder, "e_xc_solv_implicit.db"))
    name = str(Path(data[idx]["run_structure"]).stem).replace(data[idx]["metal"], "M")
    try:
        # Check if e_xc is present in the database, if not catalyst is discarded
        _ = e_xc_db.select(name=name)
    except AttributeError:
        return False, None
    converged = e_mxc(data[idx], vasp_parameters)
    if converged:
        return True, data[idx]
    else:
        return False, None


if __name__ == "__main__":
    main()
