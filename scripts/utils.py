import json
import logging
import os
import re
from itertools import combinations
from pathlib import Path

import numpy as np  # type: ignore
from ase.build import add_adsorbate # type: ignore
from ase.calculators.vasp import Vasp  # type: ignore
from ase.db import connect  # type: ignore
from ase.io import read, write  # type: ignore


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
        1 if calculation didn't stop because reached the max number of ionic steps, 0 otherwise.
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
    directory: os.PathLike, vasp_parameters: dict, calculation: str) -> bool:
    """Run the VASP calculations with NO solvation and dipole coorections.

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
    atoms = read(os.path.join(directory, "init.POSCAR"))
    mag = magmons()
    for atom in atoms:
        if atom.symbol in mag.keys():
            atom.magmom = mag[atom.symbol]
    # Begin with static calculation without dipole correction
    vasp_parameters["nelmdl"] = -8
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
        for atom in atoms:
            if atom.symbol in mag.keys():
                atom.magmom = mag[atom.symbol]
        vasp_parameters["ldipol"] = True
        vasp_parameters["idipol"] = 3
        vasp_parameters["dipol"] = atoms.get_center_of_mass(scaled=True)
        vasp_parameters["nsw"] = 999
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
                    "POTCAR",
                    "CONTCAR",
                    "OUTCAR",
                ]:
                    os.system(f"cp {f} {f}.{ext}")
                # Clean up
                for f in [
                    "CHG",
                    "CHGCAR",
                    "WAVECAR",
                    "CONTCAR",
                    "DOSCAR",
                    "EIGENVAL",
                    "PROCAR",
                ]:
                    os.remove(f)
                    converged = True
    # clean up
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
        os.system("rm %s" % f)
    if converged:
        return True
    else:
        print(f"{calculation} calculation did not converge.")
        return False
    
def operating_stability_run_vasp(
    directory: os.PathLike, vasp_parameters: dict, calculation: str, solvation: bool | None = False) -> bool:
    """Run the VASP calculations.

    Args:
        directory (Path): Path to the directory where to find the initial geometry.
        vasp_parameters (dict): Dictionary containing the VASP parameters.
        calculation (str): Type of calculation to run.
        solvation (bool): Whether to run the calculation with implicit solvation and dipole correction or not.

    Returns:
        bool: True if the calculation converged, False otherwise.
    """
    converged = False
    control_ion = 0
    os.chdir(directory)
    atoms = read(os.path.join(directory, "init.POSCAR"))
    mag = magmons()
    for atom in atoms:
        if atom.symbol in mag.keys():
            atom.magmom = mag[atom.symbol]
    control_ion = 0
    # static calculation always at the beginning
    vasp_parameters["istart"] = 0  # strart from being from scratch
    vasp_parameters["icharg"] = 2  # take superposition of atomic charge density
    vasp_parameters["nsw"] = 0
    vasp_parameters["lcharg"] = True
    vasp_parameters["lwave"] = False
    vasp_parameters["isif"] = 0
    paramscopy = vasp_parameters.copy()
    calc1 = Vasp(**paramscopy)
    atoms.set_calculator(calc1)
    atoms.get_potential_energy()
    """Check if a vasp calculation (electronic self consistance) is converged"""
    nelm = calc1.int_params["nelm"]
    control_electronic = check_electronic(nelm)
    if control_electronic == 0:
        print("Error: electronic scf")
        control_ion = 1
    else:
        ext = "preSP"
        for f in [
            "INCAR",
            "POSCAR",
            "POTCAR",
            "CONTCAR",
            "OUTCAR",
            "OSZICAR",
            "vasp.out",
        ]:
            os.system("cp %s %s.%s" % (f, f, ext))
    # relax and dipole correction calculation
    atoms = read("CONTCAR")
    for atom in atoms:
        if atom.symbol in mag.keys():
            atom.magmom = mag[atom.symbol]
    vasp_parameters["istart"] = 0  # not to read WAVECAR
    vasp_parameters["icharg"] = 1  # restrat from CHGCAR
    vasp_parameters["ldipol"] = True
    vasp_parameters["idipol"] = 3
    vasp_parameters["dipol"] = atoms.get_center_of_mass(scaled=True)
    if solvation == "implicit":
        vasp_parameters["isif"] = 0
        vasp_parameters["lwave"] = True
        vasp_parameters["lcharg"] = True
        vasp_parameters["nsw"] = 50
    else:
        vasp_parameters["isif"] = 0
        vasp_parameters["lwave"] = False
        vasp_parameters["lcharg"] = True
        vasp_parameters["nsw"] = 50
    paramscopy = vasp_parameters.copy()
    calc2 = Vasp(**paramscopy)
    atoms.set_calculator(calc2)
    atoms.get_potential_energy()
    """Check if a vasp calculation is converged"""
    nelm = calc2.int_params["nelm"]
    control_electronic = check_electronic(nelm)
    if control_electronic == 0:
        print("Error: electronic scf")
        control_ion = 1
    else:
        """Check if force converged"""
        nsw = calc2.int_params["nsw"]
        control_ion = check_ion(nsw)
        if control_ion == 1:
            if solvation == "implicit":
                ext = "preRDip"
            else:
                ext = "RDip"
            for f in [
                "INCAR",
                "POSCAR",
                "POTCAR",
                "CONTCAR",
                "OUTCAR",
                "OSZICAR",
                "vasp.out",
            ]:
                os.system("cp %s %s.%s" % (f, f, ext))
    # clean up
    for f in ["DOSCAR", "EIGENVAL", "PROCAR", "vasprun.xml"]:
        os.system("rm %s" % f)
    control_ion = 0
    os.system("cp CONTCAR.preRDip CONTCAR")
    os.system("cp WAVECAR.preRDip WAVECAR")
    # relaxation with solvation and dipole correction
    atoms = read("CONTCAR")
    for atom in atoms:
        if atom.symbol in mag.keys():
            atom.magmom = mag[atom.symbol]
    vasp_parameters["istart"] = 1  # start from WAVECAR in vaccum
    vasp_parameters["nsw"] = 500
    vasp_parameters["ldipol"] = True
    vasp_parameters["idipol"] = 3
    vasp_parameters["dipol"] = atoms.get_center_of_mass(scaled=True)
    vasp_parameters["isif"] = 0
    vasp_parameters["lwave"] = True
    vasp_parameters["lcharg"] = True
    vasp_parameters["lsol"] = True
    vasp_parameters["eb_k"] = 80
    vasp_parameters["isif"] = 0
    paramscopy = vasp_parameters.copy()
    calc4 = Vasp(**paramscopy)
    atoms.set_calculator(calc4)
    atoms.get_potential_energy()
    """Check if a vasp electronic calculation is converged"""
    nelm = calc4.int_params["nelm"]
    control_electronic = check_electronic(nelm)  # check electronic scf
    if control_electronic == 0:
        print("Error: electronic scf")
        control_ion = 1
    else:
        """Check if force converged"""
        nsw = calc4.int_params["nsw"]
        control_ion = check_ion(nsw)
        if control_ion == 1:
            ext = "RDip"
            for f in [
                "INCAR",
                "POSCAR",
                "POTCAR",
                "CONTCAR",
                "OUTCAR",
                "OSZICAR",
                "vasp.out",
            ]:
                os.system("cp %s %s.%s" % (f, f, ext))
                converged = True
    # clean up
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
        os.system("rm %s" % f)
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
       workflow_data: Dictionary containing the run structure, base directory and metal.
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
    structure = read(f'{outcar}@-1')
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

def get_xch_structs(data: list) -> None:
    """Get xch structures for running calculation.

    Args:
        data (dict): A dictionary containing the run parameters.

    Returns:
        None.
    """
    new_data = {}
    overall_idx = 0
    run_struc = data["run_structure"]
    name = Path(run_struc).stem
    metal = data["metal"]
    base_dir = data["base_dir"]
    skimmed_data = data.copy()
    del skimmed_data["metal"]
    structure = read(os.path.join(run_struc, "init.POSCAR"))
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
                                "e_xch_implicit_solv",
                                f"{name.replace(metal, metal_str)}",
                                f"{h}H_config{idx+1}",
                            )
                        )
                        run_dir.mkdir(exist_ok=True, parents=True)
                        write(f"{run_dir}/init.POSCAR", no_metal_copy, format="vasp")
                        merged_dict = {**skimmed_data.copy(), "run_dir": str(run_dir)}
                        new_data[overall_idx] = merged_dict
                        overall_idx += 1
    return new_data

def magmons():
    """All ferromagnetic or anti-ferromagnetic elements 
    and their magnetic moments multiplied
    by 1.2 already for VASP calculations.
    """
    elements = {
        "Co": 1.2, # multiplied by 1.2
        "Ti": 1.8, # multiplied by 1.2
        "Ni": 1.8, # multiplied by 1.2
        "Mn": 6.0, # multiplied by 1.2 ( Mn2+ ? can be paramgentic or anti-ferromagnetic)
        "Fe": 3.6, # multiplied by 1.2
        "V": 3.8, # multiplied by 1.2
        "Zr": 1.8, # multiplied by 1.2
        "Nb": 3.2, # multiplied by 1.2
        "Mo": 4.6, # multiplied by 1.2
        "Ta": 3.4, # multiplied by 1.2
        "W": 6.2, # multiplied by 1.2
        "Sc": 2.8, # multiplied by 1.2
        "Pt": 0.5,  # NOT multiplied by 1.2 as paramagnetic
        "Pd": 0.5, # NOT multiplied by 1.2 as paramagnetic
        "Cr": 3.9  # multiplied by 1.2
        }

    return elements