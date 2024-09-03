"""Code taken and refactored from Tiaporn"""

import os
from pathlib import Path

from ase.calculators.vasp import Vasp  # type: ignore
from typing import Optional # type: ignore
from ase.io import read  # type: ignore
from utils import (  # type: ignore
    check_electronic,
    check_ion,
    read_and_write_database,
    run_logger,
    magmons,
    check_database
)
from vasp_input import vasp_input  # type: ignore


def check_geometry(cwd: os.PathLike) -> int:
    """Check if carbon is too far from the metal site. If so, the structure is not reasonable.

    Args:
       cwd (os.PathLike): Current working directory.

    Returns:
       control (int): 0 if the structure is not reasonable, 1 otherwise.
    """
    control = 1
    folders = str(cwd).split("/")
    for n in range(len(folders)):
        if (
            folders[n] == "PBE+D3"
            or folders[n] == "RPBE+D3"
            or folders[n] == "BEEF-vdW"
        ):
            break
    metal = folders[n + 2]
    carbons = []
    atoms = read("CONTCAR")

    # check the metal site
    for atom in atoms:
        if atom.symbol == "C":
            carbons.append(atom.z)
        if atom.symbol == metal:
            metal_z = atom.z
    c_avg = sum(carbons) / len(carbons)
    difference = abs(metal_z - c_avg)
    if difference > 2.0:
        print("Error: metal-C distance", flush=True)
        control = 0
    return control


def relax_pristine(cwd: os.PathLike, data: dict) -> bool:
    """Relax the pristine structure in vacuum for adsorbate calculations further on.

    Args:
       cwd (os.PathLike): Current working directory.
       data (dict): Dictionary containing the run parameters.

    Returns:
       relaxed (bool): True if the calculation is done, False otherwise.
    """
    folders = str(cwd).split("/")
    for n in range(len(folders)):
        if (
            folders[n] == "PBE+D3"
            or folders[n] == "RPBE+D3"
            or folders[n] == "BEEF-vdW"
        ):
            break
    solvation = folders[n + 12]
    params = vasp_input()
    paramscopy = params.copy()
    calc = Vasp(**paramscopy)
    # start the vasp calculations
    err = 0
    control_ion = 0
    converged = False
    atoms = read("init.POSCAR")
    mag = magmons()
    for atom in atoms:
        if atom.symbol in mag.keys():
            atom.magmom = mag[atom.symbol]
    # static calculation always at the beginning
    params["istart"] = (
        0  # strart from being from scratch (ISTART determins whether or not to read the WAVECAR file)
    )
    params["icharg"] = 2  # take superposition of atomic charge density
    params["nsw"] = 0
    params["lcharg"] = True
    params["lwave"] = False
    params["isif"] = 0
    calc1 = calc
    paramscopy = params.copy()
    calc1 = Vasp(**paramscopy)
    atoms.set_calculator(calc1)
    atoms.get_potential_energy()
    """Check if a vasp calculation (electronic self consistance) is converged"""
    nelm = calc1.int_params["nelm"]
    control_electronic = check_electronic(nelm)
    if control_electronic == 0:
        print("Error: electronic scf")
        control_ion = 1
        err = 1
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
            "WAVECAR",
        ]:
            os.system("cp %s %s.%s" % (f, f, ext))
        # relax and dipole correction without solvation
        atoms = read("CONTCAR")
        for atom in atoms:
            if atom.symbol in mag.keys():
                atom.magmom = mag[atom.symbol]
        params["istart"] = 0
        params["icharg"] = 1
        params["ldipol"] = True
        params["idipol"] = 3
        params["dipol"] = atoms.get_center_of_mass(scaled=True)
        if solvation == "vac":
            params["isif"] = 4
            params["lwave"] = False
            params["lcharg"] = True
            params["nsw"] = 50
            params["ldipol"] = False
        elif solvation == "implicit":
            params["isif"] = 4
            params["lwave"] = True
            params["lcharg"] = True
            params["nsw"] = 50
        calc2 = calc
        paramscopy = params.copy()
        calc2 = Vasp(**paramscopy)
        atoms.set_calculator(calc2)
        atoms.get_potential_energy()
        """Check if a vasp calculation (eletronic self consistance) is converged"""
        nelm = calc2.int_params["nelm"]
        control_electronic = check_electronic(nelm)
        if control_electronic == 0:
            print("Error: electronic scf")
            control_ion = 1
            err = 1
        else:
            """Check if a reasonable geometry"""
            control_geometry = check_geometry(cwd)
            if control_geometry == 0:
                print("Error: structure disintegrates")
                control_ion = 1
                err = 1
            else:
                ext = "preRDip"
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
    if err == 0:
        control_ion = 0
        atoms = read("CONTCAR.preRDip")
        for atom in atoms:
            if atom.symbol in mag.keys():
                atom.magmom = mag[atom.symbol]
        # static calculation always at the beginning
        params["nsw"] = 0
        params["isif"] = 0
        params["lcharg"] = True
        params["ldipol"] = False
        if solvation == "implicit":
            os.system("cp WAVECAR.preRDip WAVECAR")
            params["istart"] = 1  # start from WAVECAR
            params["lwave"] = True
            params["lsol"] = True
            params["eb_k"] = 80
        else:
            params["istart"] = (
                0  # start from being from scratch (ISTART determins whether or not to read the WAVECAR file)
            )
            params["lwave"] = False
        calc3 = calc
        paramscopy = params.copy()
        calc3 = Vasp(**paramscopy)
        atoms.set_calculator(calc3)
        atoms.get_potential_energy()
        """Check if a vasp electronic calculation is converged"""
        nelm = calc3.int_params["nelm"]
        control_electronic = check_electronic(nelm)  # check electronic scf
        if control_electronic == 0:
            print("Error: electronic scf")
            control_ion = 1
            err = 1
        else:
            ext = "SP"
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
            # relax and dipole correction with solvation
            atoms = read("CONTCAR")
            for atom in atoms:
                if atom.symbol in mag.keys():
                    atom.magmom = mag[atom.symbol]
            params["icharg"] = 1
            params["nsw"] = 9999
            params["ldipol"] = True
            params["idipol"] = 3
            params["dipol"] = atoms.get_center_of_mass(scaled=True)
            params["isif"] = 0
            if solvation == "implicit":
                params["istart"] = 1
                params["lsol"] = True
                params["eb_k"] = 80
                params["lwave"] = True
                params["lcharg"] = True
            else:
                params["istart"] = 0
                params["lwave"] = False
                params["lcharg"] = True
            calc4 = calc
            paramscopy = params.copy()
            calc4 = Vasp(**paramscopy)
            atoms.set_calculator(calc4)
            atoms.get_potential_energy()
            """Check if a vasp electronic calculation is converged"""
            nelm = calc4.int_params["nelm"]
            control_electronic = check_electronic(nelm)  # check electronic scf
            if control_electronic == 0:
                print("Error: electronic scf")
                control_ion = 1
                err = 1
            else:
                """Check if a vasp ion calculation is converged"""
                control_geometry = check_geometry(
                    cwd
                )  # check matal-N and/or dopant-nb distance from CONTCAR
                if control_geometry == 0:
                    print("Error: structure disintegrates")
                    control_ion = 1
                    err = 1
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
                        outcar = Path(cwd) / "OUTCAR.RDip"
                        if "implicit" in str(cwd):
                            read_and_write_database(
                                outcar, "pristine_implicit", data
                            )
                        elif "vac" in str(cwd):
                            read_and_write_database(outcar, "pristine_vac", data)
            # clean up
            for f in [
                "CHG",
                "CHGCAR",
                "WAVECAR",
                "CONTCAR",
                "DOSCAR",
                "EIGENVAL",
                "PROCAR",
            ]:
                os.system("rm %s" % f)

    if converged:
        print("Relaxation is done")
        return True
    else:
        run_logger(
            f"Adsorbate relaxation calculation did not converge in {cwd}.", str(__file__), "error"
        )
        print(f"Adsorbate relaxation calculation is not done in {cwd}")
        raise ValueError(f"Adsorbate relaxation calculation is not done in {cwd}")


def main(**data: dict) -> tuple[bool, Optional[dict]]:
    """Run vac and implicit prisitine structures.

    Args:
       data (dict): Dictionary containing the run parameters.

    Returns:
       Perqueue return tuple.
    """
    original_dir = os.getcwd()
    prisitine_dir = Path(str(data["pristine"]))
    vac_dir = prisitine_dir / "vac" / "vasp_rx"
    implicit_dir = prisitine_dir / "implicit" / "vasp_rx"
    copy_data = data.copy()
    controls = []
    del copy_data["pq_index"]
    data["name"] = str(Path(data["run_structure"]).stem)
    name = str(Path(data["run_structure"]).stem)
    for directory in [vac_dir, implicit_dir]:
        os.chdir(directory)
        db_name = "pristine_implicit"
        print(f"Running {directory}", flush=True)
        cwd = os.getcwd()
        outcar = Path(cwd) / "OUTCAR.RDip"
        if check_database(db_name, copy_data, master=True):
            print('In master database already', flush=True)
            controls.append(True)
            break
        if os.path.exists("OUTCAR.RDip"):
            print(f"OUTCAR.RDip exists in {cwd}")
            if not check_database(db_name, copy_data, master=False):
                print(f'Writing to local database for {name}', flush=True)
                if "implicit" in str(cwd):
                    read_and_write_database(
                        outcar, "pristine_implicit", copy_data
                    )
                elif "vac" in str(cwd):
                    read_and_write_database(outcar, "pristine_vac", copy_data)
            controls.append(True)
        else:
            controls.append(relax_pristine(directory, copy_data))
    os.chdir(original_dir)
    del copy_data["pristine"]
    if all(controls):
        return True, copy_data
    else:
        return False, None


if __name__ == "__main__":
    main()
