import os
from itertools import combinations
from pathlib import Path

import numpy as np  # type: ignore
from ase.build import add_adsorbate  # type: ignore
from ase.io import read, write  # type: ignore
from perqueue.constants import INDEX_KW  # type: ignore
from utils import (  # type: ignore
    read_and_write_database,
    run_logger,
    operating_stability_run_vasp,
)
from vasp_input import vasp_input  # type: ignore
                        
def e_xch(data: dict, vasp_parameters: dict) -> bool:
    """Generate and run the input files for the e_xch calculations with solvent and dipole corrections.
    
    Args:
        data (dict): Dictionary containing the run parameters.
        
    Returns:
        True if the calculation converged, False otherwise.
    """
    run_dir = Path(data["run_dir"])
    config = Path(data["run_dir"]).stem
    structure = Path(data["run_dir"]).parent.stem
    if os.path.exists(run_dir / "OUTCAR.RDip"):
        outcar = Path(run_dir) / "OUTCAR.RDip"
        return True
    converged = operating_stability_run_vasp(
        run_dir, vasp_parameters, "e_xch_solv_implicit",
    )
    if converged:
        outcar = Path(run_dir) / "OUTCAR.RDip"
        read_and_write_database(outcar, "e_xch_solv_implicit", data)
        return True
    else:
        run_logger(
            f"e_xch_solv_implicit calculation did not converge for {config} in {structure}.",
            str(__file__),
            "error",
        )
        print(f"e_xch_solv_implicit calculation did not converge for {config} in {structure}.")
        raise ValueError(f"e_xch_solv_implicit calculation did not converge for {config} in {structure}.")


def main(**data: dict) -> tuple[bool, None]:
    """Run the hydrogen in place of metal calculations.

    Args:
        data (dict): Dictionary containing run parameters.
    Returns:
        Perqueue return tuple.
    """
    cwd = os.getcwd()
    vasp_parameters = vasp_input()
    idx, *_ = data[INDEX_KW]
    idx = str(idx)
    converged = e_xch(data[idx], vasp_parameters)
    os.chdir(cwd)
    return converged, None

if __name__ == "__main__":
    main()