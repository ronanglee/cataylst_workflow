import os
from pathlib import Path

from perqueue.constants import INDEX_KW  # type: ignore
from utils import (  # type: ignore
    check_database,
    check_for_duplicates_sql,
    operating_stability_run_vasp,
    read_and_write_database,
    run_logger,
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
    data["name"] = f"{structure}-{config}"
    print("entering e_xch", flush=True)
    master_database_dir = "/home/energy/rogle/asm_orr_rxn/master_databases"
    if check_database("e_xch_solv_implicit_without_corrections", data, master=True):
        print(
            "In master e_xch_solv_implicit_without_corrections database already",
            flush=True,
        )
        if check_for_duplicates_sql(
            f"{master_database_dir}/e_xch_vib_solv_implicit_corrections_master", data
        ):
            print(
                "In master e_xch_vib_solv_implicit_corrections database already",
                flush=True,
            )
            return True
    converged = operating_stability_run_vasp(
        run_dir, vasp_parameters, "e_xch_solv_implicit", data
    )
    if converged:
        outcar = Path(run_dir) / "OUTCAR.RDip"
        read_and_write_database(outcar, "e_xch_solv_implicit_without_corrections", data)
        return True
    else:
        run_logger(
            f"e_xch_solv_implicit calculation did not converge for {config} in {structure}.",
            str(__file__),
            "error",
        )
        print(
            f"e_xch_solv_implicit calculation did not converge for {config} in {structure}.",
            flush=True,
        )
        raise ValueError(
            f"e_xch_solv_implicit calculation did not converge for {config} in {structure}."
        )


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
