"""Code taken and refactored from Tiaporn"""

import os
from pathlib import Path

from ase.calculators.vasp import Vasp  # type: ignore
from ase.io import read, write  # type: ignore
from ase.vibrations import Vibrations  # type: ignore
from utils import (  # type: ignore
    check_electronic,
    check_for_duplicates_sql,
    check_ion,
    get_adsorbateid,
    get_vibrational_correction,
    insert_data,
    magmons,
    run_logger,
)
from vasp_input import vasp_input  # type: ignore


def adsorbsite(slab: dict, metal: str, name: str) -> tuple:
    """Find the adsorption site.

    Args:
       slab (dict): Dictionary containing the slab.
       metal (str): Metal symbol.
       name (str): Adsorption site name.

    Returns:
       tuple: x, y, z coordinates of the adsorption site."""
    if name[0] == "o" and name[6] == "m":
        print("adsorption site name = %s" % name, flush=True)
        for atom in slab:
            if atom.symbol == metal:
                metal_x = atom.x
                metal_y = atom.y
                metal_z = atom.z
        x = metal_x
        y = metal_y
        z = metal_z
    print("adsorption site (x, y, z): %.3f %.3f %.3f" % (x, y, z), flush=True)
    return (x, y, z)


def calc_vibration(cwd: os.PathLike, data: dict) -> bool:
    """Calculate the vibration of the adsorbate.

    Args:
       cwd (os.PathLike): Current working directory.
       data (dict): Dictionary containing the run parameters.

    Returns:
       return (bool): True if the calculation is successful.
    """
    folders = str(cwd).split("/")
    for n in range(len(folders)):
        if (
            folders[n] == "PBE+D3"
            or folders[n] == "RPBE+D3"
            or folders[n] == "BEEF-vdW"
        ):
            break
    metal = folders[n + 2]
    solvation = folders[n + 12]
    params = vasp_input()
    paramscopy = params.copy()
    struc_name = str(Path(data["run_structure"]).stem)
    data_base_folder = Path(data["base_dir"]) / "runs" / "databases"
    database: dict = {}
    database[struc_name] = {}
    database[struc_name]["metal"] = metal
    database[struc_name]["ads1"] = "non"
    database[struc_name]["ads2"] = "OOH"
    database["name"] = f"{list(database.keys())[0]}"
    if os.path.exists("vibration.txt"):
        if not check_for_duplicates_sql(
            os.path.join(data_base_folder, "e_ads_vib_corrections"), database
        ):
            correction = get_vibrational_correction()
            database[struc_name]["correction"] = correction
            print("Inserting into vibrational database", flush=True)
            insert_data(
                os.path.join(data_base_folder, "e_ads_vib_corrections"),
                [list(database.keys())[0], list(database.values())[0]],
            )
        return True
    err = 0
    mag = magmons()
    calc = Vasp(**paramscopy)
    if not os.path.exists(os.path.join(cwd, "OUTCAR.RDip")):
        init_atoms = read(Path(cwd).parent / "vasp_rx" / "CONTCAR.RDip")
        write("initial_ads.POSCAR", init_atoms)
        atoms = read("initial_ads.POSCAR")
        mag = magmons()
        for atom in atoms:
            if atom.symbol in mag.keys():
                atom.magmom = mag[atom.symbol]
        # static calculation always at the beginning
        params["istart"] = 0  # strart from being from scratch
        params["icharg"] = 2  # take superposition of atomic charge density
        params["nsw"] = 0
        params["lcharg"] = True
        params["lwave"] = False
        params["isif"] = 0
        calc1 = calc
        paramscopy = params.copy()
        calc1 = Vasp(**paramscopy)
        atoms.set_calculator(calc1)
        print("static calculation", flush=True)
        atoms.get_potential_energy()
        """Check if a vasp calculation (electronic self consistance) is converged"""
        nelm = calc1.int_params["nelm"]
        control_electronic = check_electronic(nelm)
        if control_electronic == 0:
            print("Error: electronic scf", flush=True)
            err = 1
        else:
            adsorbate_id = get_adsorbateid("initial_ads.POSCAR", len("OOH"))
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
            # dipole correction calculation
            atoms = read("CONTCAR")
            for atom in atoms:
                if atom.symbol in mag.keys():
                    atom.magmom = mag[atom.symbol]
            params["istart"] = 0  # not to read WAVECAR
            params["icharg"] = 1  # restrat from CHGCAR
            params["ldipol"] = True
            params["idipol"] = 3
            params["dipol"] = atoms.get_center_of_mass(scaled=True)
            params["isif"] = 0
            params["lwave"] = True
            params["lcharg"] = True
            params["nsw"] = 0
            calc2 = calc
            paramscopy = params.copy()
            calc2 = Vasp(**paramscopy)
            atoms.set_calculator(calc2)
            print("dipole calculation", flush=True)
            atoms.get_potential_energy()
            """Check if a vasp calculation is converged"""
            nelm = calc2.int_params["nelm"]
            control_electronic = check_electronic(nelm)  # check electronic scf
            if control_electronic == 0:
                print("Error: electronic scf", flush=True)
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
                atoms = read("CONTCAR.preRDip")
                for atom in atoms:
                    if atom.symbol in mag.keys():
                        atom.magmom = mag[atom.symbol]
                # dipole calculation + solvation
                os.system("cp WAVECAR.preRDip WAVECAR")
                params["istart"] = 1
                params["nsw"] = 9999
                params["ldipol"] = True
                params["idipol"] = 3
                params["dipol"] = atoms.get_center_of_mass(scaled=True)
                params["isif"] = 0
                params["lwave"] = True
                params["lcharg"] = True
                params["lsol"] = True
                params["eb_k"] = 80
                calc4 = calc
                paramscopy = params.copy()
                calc4 = Vasp(**paramscopy)
                atoms.set_calculator(calc4)
                print("dipole correction calculation", flush=True)
                atoms.get_potential_energy()
                """Check if a vasp electronic calculation is converged"""
                nelm = calc4.int_params["nelm"]
                control_electronic = check_electronic(nelm)  # check electronic scf
                if control_electronic == 0:
                    print("Error: electronic scf", flush=True)
                    err = 1
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
                        "WAVECAR",
                    ]:
                        os.system("cp %s %s.%s" % (f, f, ext))
                    vib_done = True
    # vibration calculation
    if os.path.exists("OUTCAR.RDip") or vib_done:
        if err == 0:
            print("#vibration calculation", flush=True)
            atoms = read("CONTCAR.RDip")
            for atom in atoms:
                if atom.symbol in mag.keys():
                    atom.magmom = mag[atom.symbol]
            os.system("cp WAVECAR.RDip WAVECAR")
            params["istart"] = 1
            params["ldipol"] = True
            params["idipol"] = 3
            params["dipol"] = atoms.get_center_of_mass(scaled=True)
            params["isif"] = 0
            params["nsw"] = 999
            if solvation == "implicit":
                params["lsol"] = True
                params["eb_k"] = 80
            else:
                params["lsol"] = False
            calc7 = calc
            paramscopy = params.copy()
            calc7 = Vasp(**paramscopy)
            atoms.set_calculator(calc7)
            atoms.get_potential_energy()
            """Check if a vasp electronic calculation is converged"""
            nelm = calc7.int_params["nelm"]
            control_electronic = check_electronic(nelm)  # check electronic scf
            if control_electronic == 0:
                print("Error: electronic scf", flush=True)
                err = 1
            else:
                init_poscar = read(os.path.join(data["run_structure"], "init.POSCAR"))
                indices = [atom.index for atom in init_poscar]
                current_indices = [atom.index for atom in atoms]
                len_vib = len(current_indices) - len(indices)
                adsorbate_id = get_adsorbateid("initial_ads.POSCAR", len_vib)
                nsw = calc7.int_params["nsw"]
                control_ion = check_ion(nsw)
                if control_ion == 1:
                    ext = "vib"
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
                    vibindices = adsorbate_id
                    vib = Vibrations(
                        atoms, indices=vibindices, name="vib", delta=0.01, nfree=2
                    )
                    vib.run()
                    vib.get_energies()
                    vib.summary(log="vibration.txt")
                    for i in range(3 * len(vibindices)):
                        vib.write_mode(i)
                    correction = get_vibrational_correction()
                    database[struc_name]["correction"] = correction
                    insert_data(
                        os.path.join(data_base_folder, "e_ads_vib_corrections"),
                        [list(database.keys())[0], list(database.values())[0]],
                    )
                    # clean up
                    for f in [
                        "CHG",
                        "CHGCAR",
                        "WAVECAR",
                        "DOSCAR",
                        "EIGENVAL",
                        "PROCAR",
                    ]:
                        os.system("rm %s" % f)
                    converged = True

    if converged:
        return True
    else:
        run_logger(
            f"Vibration calculation did not converge in {cwd}.", str(__file__), "error"
        )
        print(f"Vibration calculation did not converge in {cwd}.", flush=True)
        raise ValueError(f"Vibration calculation did not converge in {cwd}.")


def main(**data: dict) -> tuple[bool, dict]:
    """Run the vibration part of the workflow.

    Args:
       data (dict): Dictionary containing the run parameters.

    Returns:
       Perqueue return tuple.
    """
    cwd = os.getcwd()
    vib_dir = Path(str(data["adsorbate"])) / "implicit" / "vibration"
    os.chdir(vib_dir)
    master_database_dir = "/home/energy/rogle/asm_orr_rxn/master_databases"
    if check_for_duplicates_sql(
        f"{master_database_dir}/e_ads_vib_corrections_master", data
    ):
        print("In master database already", flush=True)
        control_vibration = True
    else:
        control_vibration = calc_vibration(vib_dir, data)
    os.chdir(cwd)
    del data["adsorbate"]
    return control_vibration, data


if __name__ == "__main__":
    main()
