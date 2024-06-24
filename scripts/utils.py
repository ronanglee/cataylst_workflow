import json
import logging
import os
import re
from pathlib import Path

import numpy as np  # type: ignore
from ase.calculators.vasp import Vasp  # type: ignore
from ase.db import connect  # type: ignore
from ase.io import read  # type: ignore


def check_electronic(set_nelm: int) -> int:
    """Open OUTCAR and check NELM

    Args:
       set_nelm (int): NELM in INCAR

    Returns:
       1 if calculation didn't stop because reached the max number of electronic steps, 0 otherwise.
    """
    regexp = re.compile("Iteration\s+(?P<nsw>\d*)\(\s+(?P<nelm>\d*)\)")  # noqa
    max_nelm = 0
    with open("OUTCAR") as f:
        for line in f:
            if "Iteration" in line:
                m = regexp.search(line)  # type: ignore
                max_nelm = int(m.groupdict().get("nelm"))  # type: ignore
    if set_nelm > max_nelm:
        return 1
    else:
        return 0


def check_ion(set_nsw: int) -> int:
    """Open OUTCAR and check NSW

    Args:
       set_nsw (int): NSW in INCAR

    Returns:
       control (int): 1 if calculation didn't stop because reached the max number of ionic steps, 0 otherwise.
    """
    regexp = re.compile("Iteration\s+(?P<nsw>\d*)\(\s+(?P<nelm>\d*)\)")  # noqa
    max_nsw = 0
    with open("OUTCAR") as f:
        for line in f:
            if "Iteration" in line:
                m = regexp.search(line)  # type: ignore
                max_nsw = int(m.groupdict().get("nsw"))  # type: ignore
    if set_nsw > max_nsw:
        return 1
    else:
        return 0


def synthesis_stability_run_vasp(
    directory: os.PathLike, vasp_parameters: dict, calculation: str
) -> bool:
    """Run the VASP calculations.

    Args:
        directory (Path): Path to the directory where to find the initial geometry.
        vasp_parameters (dict): Dictionary containing the VASP parameters.
        calculation (str): Type of calculation to run.

    Returns:
        bool: True if the calculation converged, False otherwise.
    """
    converged = False
    control_ion = 0
    os.chdir(directory)
    vasp_parameters["amin"] = 0.1
    vasp_parameters["amix"] = 0.02
    vasp_parameters["bmix"] = 1.0
    if vasp_parameters["ispin"] == 2:  # Check if spin polarized
        vasp_parameters["amix_mag"] = 0.08
        vasp_parameters["bmix_mag"] = 1.0
    atoms = read(os.path.join(directory, "init.POSCAR"))
    # Begin with static calculation without dipole correction
    vasp_parameters["nelmdl"] = -8
    vasp_parameters["ibrion"] = 2
    vasp_parameters["nsw"] = 0
    paramscopy = vasp_parameters.copy()
    calc = Vasp(**paramscopy)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    nelm = calc.int_params["nelm"]
    control_electronic = check_electronic(nelm)
    if control_electronic == 0:
        print("Error: electronic scf")
        control_ion = 1
    else:
        """Check if force converged"""
        nsw = calc.int_params["nsw"]
        control_ion = check_ion(nsw)
        if control_ion == 1:
            ext = "SP"
            for f in [
                "OSZICAR",
                "vasp.out",
                "INCAR",
                "POSCAR",
                "POSCAR",
                "POTCAR",
                "CONTCAR",
                "OUTCAR",
            ]:
                os.system(f"cp {f} {f}.{ext}")
            # Restart Calculation from CHGCAR include dipole correction
            vasp_parameters["icharg"] = 1
            vasp_parameters["ldipol"] = True
            vasp_parameters["idipol"] = 3
            vasp_parameters["dipol"] = atoms.get_center_of_mass(scaled=True)
            vasp_parameters["nsw"] = 0
            vasp_parameters["lwave"] = True
            paramscopy = vasp_parameters.copy()
            calc = Vasp(**paramscopy)
            atoms.set_calculator(calc)
            atoms.get_potential_energy()
            """Check if a vasp calculation (eletronic self consistance) is converged"""
            nelm = calc.int_params["nelm"]
            control_electronic = check_electronic(nelm)
            if control_electronic == 0:
                print("Error: electronic scf")
                control_ion = 1
            else:
                """Check if force converged"""
                nsw = calc.int_params["nsw"]
                control_ion = check_ion(nsw)
                if control_ion == 1:
                    ext = "opt"
                    for f in [
                        "OSZICAR",
                        "vasp.out",
                        "INCAR",
                        "POSCAR",
                        "POSCAR",
                        "POTCAR",
                        "CONTCAR",
                        "OUTCAR",
                        "WAVECAR",
                    ]:
                        os.system(f"cp {f} {f}.{ext}")
                    # Clean up
                    for f in [
                        "CHG",
                        "CHGCAR",
                        "WAVECAR",
                        "POSCAR",
                        "CONTCAR",
                        "DOSCAR",
                        "EIGENVAL",
                        "PROCAR",
                    ]:
                        os.remove(f)
                        converged = True
    if converged:
        return True
    else:
        print(f"{calculation} calculation did not converge.")
        return False


def gather_structs(data: dict) -> dict:
    """Gather all the structures for next part of workflow.

    Args:
       data (dict): Dictionary containing the run structure, base directory and metal.

    Returns:
       dict: Dictionary containing the run structure, base directory and metal.
    """
    workflow_data = {}
    run_structures = data["all_run_structures"]
    metal = data["metal"]
    dopants = data["dopants"]
    carbon_structure = data["carbon_structure"]
    wanted_structures = []
    for struc_path in run_structures:
        structure = read(Path(struc_path) / Path("init.POSCAR"))
        struc_metal: str = "Pt"
        if structure[35].symbol == struc_metal:
            cs = "armchair"
        elif structure[45].symbol == struc_metal or structure[38].symbol == struc_metal:
            cs = "zigzag"
        elif structure[62].symbol == struc_metal:
            cs = "bulk"
        if carbon_structure == cs:
            wanted_structures.append(struc_path)
    idx = 0
    for struc_path in wanted_structures:
        structure = read(Path(struc_path) / Path("init.POSCAR"))
        if dopants and "B" in structure.get_chemical_symbols():
            for dopant in dopants:
                data = {
                    "base_dir": data["base_dir"],
                    "run_structure": str(
                        struc_path.replace("B", dopant).replace("Pt", metal)
                    ),
                    "carbon_structure": carbon_structure,
                    "metal": metal,
                    "dopant": dopant,
                }
                workflow_data[idx] = data
                idx += 1
        else:
            workflow_data[idx] = {
                "base_dir": data["base_dir"],
                "run_structure": str(struc_path).replace("Pt", metal),
                "carbon_structure": carbon_structure,
                "metal": metal,
                "dopant": "",
            }
            idx += 1
    return workflow_data


def get_vibrational_correction() -> float:
    """Calculate vibrational correction.

    Returns:
        float: Vibrational correction.
    """
    kt = 0.02568
    rt = 298.15
    vib_file = open("vibration.txt")
    zpe_tot = 0
    u_tot = 0
    ts_tot = 0
    for line in vib_file:
        pattern = re.compile(r"\s*(\d+)\s+(?P<energy>\d+\.\d+)\s+(\d+\.\d+)\s*")
        match = pattern.match(line)
        if match:
            finder = pattern.search(line)
            e = float(finder.groupdict()["energy"]) / 1000  # type: ignore
            u = e / (np.exp(e / kt) - 1)
            zpe_tot += 0.5 * e  # type: ignore
            u_tot += e / (np.exp(e / kt) - 1)
            ts_tot += rt * (kt / rt) * (u / kt - np.log(1 - np.exp(-e / kt)))
    correction = zpe_tot + u_tot - ts_tot
    return correction


def load_data(file_path: str) -> dict:
    """Load data from a json file.

    Args:
        file_path (str): Path to the json file.

    Returns:
        Data from the json file.
    """
    try:
        with open(file_path, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        return {}


def save_data(file_path: str, data: dict) -> None:
    """Save data to a json file.

    Args:
        file_path (str): Path to the json file.
        data (dict): Data to be saved.

    Returns:
        None
    """
    with open(file_path, "w") as file:
        json.dump(data, file, indent=4)


def add_entry(file_path: str, new_entry: dict) -> None:
    """Add a new entry to a json file.

    Args:
        file_path (str): Path to the json file.
        new_entry (dict): New entry to be added.

    Returns:
        None
    """
    data = load_data(file_path)
    data.update(new_entry)
    save_data(file_path, data)


def read_and_write_database(outcar: os.PathLike, database: str, data: dict) -> None:
    """Read and write the database.

    Args:
       outcar (Path): Path to the OUTCAR file.
       name (str): Name of the structure.
       database (str): Name of the database.
       data (dict): Dictionary containing the data to be written to the database.

    Returns:
       None

    """
    data_base_folder = Path(data["base_dir"]) / "runs" / "databases"
    db = connect(data_base_folder / f"{database}.db")
    structure = read(outcar)
    db.write(structure, key_value_pairs={**data})


def run_logger(msg: str, filename: str, type_msg: str) -> None:
    """Run logger for tracking errors and discarded values.

    Args:
        msg (str): Error message.
        filename (str): Name of the file where the error occurred.
        type_msg (str): Type of the message.

    Returns:
        None
    """
    base_dir = Path(__file__).parent
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(base_dir / "workflow.log", mode="a"),
            logging.StreamHandler(),
        ],
    )

    logger = logging.getLogger(filename)
    if type_msg == "error":
        logger.error(msg)
    else:
        logger.info(msg)
