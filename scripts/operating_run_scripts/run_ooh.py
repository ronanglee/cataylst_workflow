"""Code taken and refactored from Tiaporn"""

import json
import os
from pathlib import Path

from ase.io import read, write  # type: ignore


def find_structurers(base_dir: os.PathLike) -> None:
    """Find the filtered structures for the operating stability calculations.

    Args:
        base_dir (os.PathLike): Path to the base directory.

    Returns:
        None
    """
    run_dir = Path(base_dir) / "old_good_run"
    with open(run_dir / "seperated_operating_stability.json", "r") as file:
        return json.load(file)


def main(**data: dict) -> tuple[bool, dict]:
    """Generate input files for OOH adsorbed to catalyst.

    Args:
        data (dict): Dictionary containing the base directory, run structure, carbon structure and metal.

    Returns:
        Perqueue return tuple.
    """
    calc = "non_and_adsorbate"
    xc = "BEEF-vdW"
    dopant = "non"
    dopantsite = "non"
    ads1 = "non"
    ads1site = "non"
    ads1orient = "non"
    ads2 = ["non", "OOH"]
    ads2site = ["non", "ontop_metal"]
    ads2orient = ["non", "1"]
    solvation = ["implicit", "vac"]
    metal = data["metal"]
    prototype = Path(str(data["run_structure"])).stem
    base_dir = data["base_dir"]

    for ad2, ad2site, ad2orient in zip(ads2, ads2site, ads2orient):
        for solvent in solvation:
            if solvent == "vac" and ad2 != "non":
                continue  # Skip the vac calculation for OOH as not needed
            if (
                solvent == "implicit" and ad2 != "non"
            ):  # Do vibration calculation for OOH
                run_folders = [
                    "vibration",
                    "vasp_rx",
                ]
            else:
                run_folders = ["vasp_rx"]
            for run_folder in run_folders:
                folder_structure = Path(
                    os.path.join(  # type: ignore
                        str(base_dir),
                        "runs",
                        "operating_stability",
                        calc,
                        xc,
                        prototype,
                        metal,
                        "1",
                        dopant,
                        dopantsite,
                        ads1,
                        ads1site,
                        ads1orient,
                        ad2,
                        ad2site,
                        ad2orient,
                        solvent,
                        run_folder,
                    )
                )  # type: ignore
                folder_structure.mkdir(parents=True, exist_ok=True)
                if ad2 == "non":
                    data["pristine"] = str(folder_structure.parent.parent)  # type: ignore
                elif ad2 != "non" and run_folder == "vasp_rx":
                    data["adsorbate"] = str(folder_structure.parent.parent)  # type: ignore
                elif ad2 != "OOH" and run_folder == "vibration":
                    data["vibration"] = str(folder_structure.parent.parent)  # type: ignore
                slab = read(os.path.join(data["run_structure"], "init.POSCAR"))  # type: ignore
                for atom in slab:
                    if atom.symbol == "Pt":
                        atom.symbol = metal
                write(
                    os.path.join(folder_structure, "init.POSCAR"), slab, format="vasp"
                )
    return True, data


if __name__ == "__main__":
    main()
