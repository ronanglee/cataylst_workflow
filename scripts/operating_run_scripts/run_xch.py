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
    synthesis_stability_run_vasp,
)
from vasp_input import vasp_input  # type: ignore


def e_xch(data: dict, vasp_parameters: dict) -> None:
    """Run the xch calculations for operating stability.

    Args:
        data (dict): A dictionary containing the run parameters.
        vasp_parameters (dict): A dictionary containing the VASP parameters.

    Returns:
        None.
    """
    run_struc = data["run_structure"]
    name = Path(run_struc).stem
    metal = data["metal"]
    base_dir = data["base_dir"]
    structure = read(os.path.join(run_struc, "init.POSCAR"))
    skimmed_data = data.copy()
    del skimmed_data["metal"]
    m_x = [atom.x for atom in structure if atom.symbol == metal]
    m_y = [atom.y for atom in structure if atom.symbol == metal]
    indices = {}
    for atom in structure:
        if atom.symbol != metal:
            atom_x = atom.x
            atom_y = atom.y
            distance = list(np.sqrt((m_x - atom_x) ** 2 + (m_y - atom_y) ** 2))[0]
            if distance < 2.2:
                indices[atom.index] = atom.symbol
    no_metal = structure.copy()
    del no_metal[[atom.index for atom in structure if atom.symbol == metal]]
    for h in range(1, len(indices) + 1):
        no_metal_copy = no_metal.copy()
        idicess = []
        for atom in no_metal:
            if atom.index in indices.keys():
                idicess.append(atom.index)
                if len(idicess) == h:
                    perm = combinations(indices.keys(), h)
                    idxx = 0
                    combos = []
                    for p in perm:
                        if idxx < 5 - h:
                            combos.append(p)
                        idxx += 1
                    for idx, combo in enumerate(combos):
                        no_metal_copy = no_metal.copy()
                        for c in combo:
                            add_adsorbate(
                                no_metal_copy,
                                "H",
                                height=1.0,
                                position=(no_metal[c].x, no_metal[c].y),
                            )
                        metal_str = "M"  # so f-string wont give error
                        run_dir = Path(
                            os.path.join(
                                base_dir,
                                "runs",
                                "operating_stability",
                                "e_xch",
                                f"{name.replace(metal, metal_str)}",
                                f"{h}H_config{idx+1}",
                            )
                        )
                        run_dir.mkdir(exist_ok=True, parents=True)
                        write(f"{run_dir}/init.POSCAR", no_metal_copy, format="vasp")
                        converged = synthesis_stability_run_vasp(
                            run_dir, vasp_parameters, "e_xch"
                        )
                        if converged:
                            outcar = Path(run_dir) / "OUTCAR.opt"
                            skimmed_data["name"] = str(
                                Path(data["run_structure"]).stem
                            ).replace(metal, f"config{idx}_")
                            read_and_write_database(outcar, "e_xch", skimmed_data)
                        else:
                            run_logger(
                                f"config{idx} calculation did not converge.",
                                str(__file__),
                                "error",
                            )
                            print(f"config{idx} calculation did not converge.")


def main(**data: dict) -> tuple[bool, None]:
    """Run the synthesis stability part of the workflow.
    g_a = e_xc + e_m_on_c - e_c - e_mxc
    g_d = e_m + e_xc - e_mxc
    where x is the dopant, c is carbon and m is the metal.

    Args:
        data (dict): Dictionary containing run parameters.
    Returns:
        Perqueue return tuple.
    """
    cwd = os.getcwd()
    vasp_parameters = vasp_input()
    idx, *_ = data[INDEX_KW]
    idx = str(idx)
    e_xch(data[idx], vasp_parameters)
    os.chdir(cwd)
    return True, None


if __name__ == "__main__":
    main()
