import math
import os
from pathlib import Path
from typing import Optional  # type: ignore

from ase.build import add_adsorbate  # type: ignore
from ase.calculators.vasp import Vasp  # type: ignore
from ase.db import connect  # type: ignore
from ase.io import read, write  # type: ignore
from utils import (  # type: ignore
    check_ase_database,
    check_electronic,
    check_ion,
    magmons,
    read_and_write_database,
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


def place_adsorbate(cwd: os.PathLike, data: dict) -> None:
    """Place the adsorbate on the slab.

    Args:
        cwd (os.PathLike): Current working directory.
        data (dict): Dictionary containing the run parameters.

    Returns:
        None
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
    ads2 = folders[n + 9]
    ads2site = folders[n + 10]
    ads2orient = folders[n + 11]
    solvation = folders[n + 12]
    pristine_vac = Path(data["base_dir"]) / "runs" / "databases"
    struc_name = str(Path(data["run_structure"]).stem)
    slab_database = connect(os.path.join(pristine_vac, "pristine_vac.db"))
    slab = slab_database.get_atoms(name=struc_name)
    ooh_db_dir = Path(data["base_dir"]) / "stock_databases"
    ooh_db = connect(os.path.join(ooh_db_dir, "ooh_ads.db"))
    solv_ads = f"{ads2}-{ads2orient}_{solvation}"
    adsorb2 = ooh_db.get_atoms(name=solv_ads)
    name = ads2site
    x2, y2, _ = adsorbsite(slab, metal, name)
    if solvation == "vac" or solvation == "implicit":
        molid0 = [
            "NC",
            "SCN",
            "O",
            "H2O",
            "H2O2",
            "COOH",
            "H",
            "F",
            "Cl",
            "Br",
            "HCN",
            "HF",
            "HCl",
            "HBr",
            "CO2",
            "CHO",
            "CO2near",
            "COH",
            "CHOH",
            "CH2OH",
            "CH2",
            "CH3",
        ]
        molid1 = [
            "CN",
            "NCS",
            "OH",
            "OOH",
            "NO2",
            "O2",
            "CO",
            "HCOO",
            "ON",
            "NO",
            "HCOOH",
            "HNO3",
            "HNO2",
            "HSCN",
            "OCH2",
            "OCH3",
        ]
        molid2 = ["OCOH", "HCO3", "CO3", "NO3", "H2CO3", "H2CO3tt"]
        molid3 = ["CO3"]
        molid4 = [
            "HSO4",
            "SO4",
            "H2PO4",
            "H2SO4",
            "HPO4",
            "PO4",
            "H3PO4",
            "ClO4",
            "HClO4",
        ]
        molid5 = ["H2CO3ct"]
        if ads2 in molid0:
            mol_index2 = 0
        elif ads2 in molid1:
            if ads2 == "NO2":
                if ads2orient == 1:
                    mol_index2 = 1
                else:
                    mol_index2 = 0
            mol_index2 = 1
        elif ads2 in molid2:
            mol_index2 = 2
        elif ads2 in molid3:
            mol_index2 = 3
        elif ads2 in molid4:
            mol_index2 = 4
        elif ads2 in molid5:
            mol_index2 = 5

        if ads2orient == "2":
            if ads2 == "2OH" or ads2 == "H2O2":
                x2 = x2 + 0.8
            elif ads2 == "O2":
                x2 = x2 - 0.65
            elif ads2 == "OOH":
                x2 = x2 + 0.65
            elif ads2 == "NO2":
                x2 = x2 - 1.0

    if ads2 == "2OH" and ads2orient == "2":
        h = 1.6
    elif ads2 == "O2" and ads2orient == "2":
        h = 1.6
    elif ads2 == "OOH" and ads2orient == "2":
        h = 1.6
    elif ads2 == "H" and ads2orient == "1":
        h = 1.6
    elif ads2 == "CO3":
        h = 4.0
    elif ads2 == "ON":
        h = 3.5
    elif ads2 == "Cl" or ads2 == "Br" or ads2 == "HCl" or ads2 == "HBr":
        h = 2.5
    elif ads2 == "CO2near":
        h = 1.8
    elif ads2 == "NO2" and ads2orient == "2":
        h = 2.6
    else:
        h = 2.1
    add_adsorbate(
        slab, adsorb2, h, position=(x2, y2), offset=None, mol_index=mol_index2
    )
    write("initial_adsorbate.POSCAR", slab, format="vasp")


def get_adsorbateid(cwd: os.PathLike) -> list:
    """Get the adsorbate id.

    Args:
        cwd (os.PathLike): Current working directory.

    Returns:
        adsorbate_id (int): Adsorbate atoms IDs.
    """
    folders = str(cwd).split("/")
    for n in range(len(folders)):
        if (
            folders[n] == "PBE+D3"
            or folders[n] == "RPBE+D3"
            or folders[n] == "BEEF-vdW"
        ):
            break
    ads2 = folders[n + 9]
    system = read("OUTCAR@-1")
    inital_distances = []
    for idx, atom in enumerate(system):
        inital_distances.append((idx, atom.z))

    sorted_list = sorted(inital_distances, key=lambda t: t[1])[-len(ads2) :]  # noqa
    adsorbate_id = [i[0] for i in sorted_list]
    return adsorbate_id


def check_geometry(cwd: os.PathLike, adsorbate_id: list) -> int:
    """Check if the geometry is reasonable.

    Args:
        cwd (os.PathLike): Current working directory.
        adsorbate_id (list): Adsorbate atoms IDs.

    Returns:
        control (int): 0 if the geometry is not reasonable, 1 otherwise.
    """
    print("Checking geometry...", flush=True)
    control = 1
    folders = str(cwd).split("/")
    for n in range(len(folders)):
        if (
            folders[n] == "PBE+D3"
            or folders[n] == "RPBE+D3"
            or folders[n] == "BEEF-vdW"
        ):
            break
    ads2 = folders[n + 9]
    metal = folders[n + 2]
    ads2site = folders[n + 10]
    name = ads2site
    system = read("OUTCAR@-1")
    initial_system = read("initial_adsorbate.POSCAR")
    x_ads, y_ads, z_ads = adsorbsite(initial_system, metal, name)
    ads_dist = []
    atomsx = []
    atomsy = []
    atomsz = []
    symbol = []
    for atom in system:
        for i in range(len(adsorbate_id)):
            if atom.index == adsorbate_id[i]:
                symbol.append(atom.symbol)
                atomsx.append(atom.x)
                atomsy.append(atom.y)
                atomsz.append(atom.z)
                dx = abs(atom.x - x_ads)
                dy = abs(atom.y - y_ads)
                dz = abs(atom.z - z_ads)
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                ads_dist.append((atom.index, dist, atom.x, atom.y, atom.z))
    min_arr = min(ads_dist, key=lambda t: t[1])
    for atom in system:
        if atom.index == min_arr[0]:
            xx = atom.x
            yy = atom.y
            zz = atom.z
    dx = abs(xx - x_ads)
    dy = abs(yy - y_ads)
    dz = zz - z_ads
    print(
        "Adsorbate (index, x, y, z): %d %.3f %.3f %.3f" % (min_arr[0], xx, yy, zz),
        flush=True,
    )
    rxy = math.sqrt(dx * dx + dy * dy)

    """ 1. check adsorbate """
    if rxy > 1.5:
        print(
            "The adsorbate shifts from adsorption site (xy > 1.5 Angstrom)", flush=True
        )
        control = 0
    elif dz > 3.5:
        print("The adsorbate deadsorps (z > 3.5 Angstrom)", flush=True)
        control = 0

    """ 2. check if the adsorbate dissociate """
    if ads2 == "OOH":
        n = symbol.count("O")
        if n == 2:
            ox = []
            oy = []
            oz = []
            for i in range(len(adsorbate_id)):
                if symbol[i] == "O":
                    ox.append(atomsx[i])
                    oy.append(atomsy[i])
                    oz.append(atomsz[i])
            oo_dist = math.sqrt(
                (ox[0] - ox[1]) ** 2 + (oy[0] - oy[0]) ** 2 + (oz[0] - oz[1]) ** 2
            )
            if oo_dist > 2.0:
                print("The adsorbate dissociates", flush=True)
                control = 0

    """ 3. check metal atom """
    carbons = []
    for atom in system:
        if atom.symbol == "C":
            carbons.append(atom.z)
        if atom.symbol == metal:
            metal_z = atom.z
    c_avg = sum(carbons) / len(carbons)
    difference = abs(metal_z - c_avg)
    if difference > 2.0:
        print("The surface reconstructs", flush=True)
        control = 0
    return control


def relax_adsorbate(cwd: os.PathLike, data: dict) -> bool:
    """Relax the adsorbate.

    Args:
        cwd (os.PathLike): Current working directory.
        data (dict): Dictionary containing the run parameters.

    Returns:
        return (bool) True if the calculation is done, False otherwise.
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
    converged = False
    qstep = 1
    if solvation == "implicit":
        qstep = 2
    print("#vasp calculation = %d" % qstep, flush=True)
    err = 0
    for ii in range(qstep):
        if ii == 0:
            control_ion = 0
            atoms = read("initial_adsorbate.POSCAR")
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
            print("(%d) static calculation" % ii, flush=True)
            atoms.get_potential_energy()
            """Check if a vasp calculation (electronic self consistance) is converged"""
            nelm = calc1.int_params["nelm"]
            control_electronic = check_electronic(nelm)
            if control_electronic == 0:
                print("Error: electronic scf", flush=True)
                control_ion = 1
                err = 1
            else:
                adsorbate_id = get_adsorbateid(cwd)
                slab = read("CONTCAR")
                for i in range(len(adsorbate_id)):
                    for atom in slab:
                        if atom.index == adsorbate_id[i]:
                            print(
                                "Adsorbate (index, x, y, z): %d %.3f %.3f %.3f"
                                % (atom.index, atom.x, atom.y, atom.z)
                            )
                if solvation == "implicit":
                    ext = "preSP"
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
            # relax and dipole correction calculation
            atoms = read("CONTCAR")
            for atom in atoms:
                if atom.symbol in mag.keys():
                    atom.magmom = mag[atom.symbol]
            params["istart"] = 0  # not to read WAVECAR
            params["icharg"] = 1  # restart from CHGCAR
            params["ldipol"] = True
            params["idipol"] = 3
            params["dipol"] = atoms.get_center_of_mass(scaled=True)
            if solvation == "implicit":
                params["isif"] = 0
                params["lwave"] = True
                params["lcharg"] = True
                params["nsw"] = 50
            else:
                params["isif"] = 0
                params["lwave"] = False
                params["lcharg"] = True
                params["nsw"] = 50
            calc2 = calc
            paramscopy = params.copy()
            calc2 = Vasp(**paramscopy)
            atoms.set_calculator(calc2)
            print("(%d) optimization calculation" % ii, flush=True)
            atoms.get_potential_energy()
            """Check if a vasp calculation is converged"""
            nelm = calc2.int_params["nelm"]
            control_electronic = check_electronic(nelm)
            if control_electronic == 0:
                print("Error: electronic scf", flush=True)
                control_ion = 1
                err = 1
            else:
                control_geometry = check_geometry(cwd, adsorbate_id)
                if control_geometry == 0:
                    print("Error: structure disintegrates", flush=True)
                    control_ion = 1
                    err = 1
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
        if ii == 1 and err == 0:
            control_ion = 0
            os.system("cp CONTCAR.preRDip CONTCAR")
            os.system("cp WAVECAR.preRDip WAVECAR")
            # relaxation with solvation and dipole correction
            atoms = read("CONTCAR")
            for atom in atoms:
                if atom.symbol in mag.keys():
                    atom.magmom = mag[atom.symbol]
            params["istart"] = 1  # start from WAVECAR in vaccum
            params["nsw"] = 500
            params["ldipol"] = True
            params["idipol"] = 3
            params["dipol"] = atoms.get_center_of_mass(scaled=True)
            params["isif"] = 0
            params["lwave"] = True
            params["lcharg"] = True
            params["lsol"] = True
            params["eb_k"] = 80
            params["isif"] = 0
            calc4 = calc
            paramscopy = params.copy()
            calc4 = Vasp(**paramscopy)
            atoms.set_calculator(calc4)
            print("(%d) optimization calculation" % ii, flush=True)
            atoms.get_potential_energy()
            """Check if a vasp electronic calculation is converged"""
            nelm = calc4.int_params["nelm"]
            control_electronic = check_electronic(nelm)  # check electronic scf
            if control_electronic == 0:
                print("Error: electronic scf", flush=True)
                control_ion = 1
                err = 1
            else:
                """Check if a vasp ion calculation is converged"""
                control_geometry = check_geometry(cwd, adsorbate_id)
                if control_geometry == 0:
                    print("Error: structure disintegrates", flush=True)
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
            outcar = Path(cwd) / "OUTCAR.RDip"
            read_and_write_database(outcar, "e_adsorbate_without_corrections", data)
            converged = True

    if converged:
        print("Relaxation is done", flush=True)
        return True
    else:
        run_logger(
            f"Adsorbate relaxation calculation is not converged in {cwd}.",
            str(__file__),
            "error",
        )
        print(
            f"Adsorbate relaxation calculation is not converged in {cwd}.", flush=True
        )
        raise ValueError(f"Adsorbate relaxation calculation is not converged in {cwd}.")


def main(**data: dict) -> tuple[bool, Optional[dict]]:
    """Run vac and implicit prisitine structures.

    Args:
        data (dict): Dictionary containing the run parameters.

    Returns:
        Perqueue return tuple.
    """
    cwd = os.getcwd()
    ooh_dir = Path(str(data["adsorbate"])) / "implicit" / "vasp_rx"
    os.chdir(ooh_dir)
    print(f"Relaxing adsorbate in {ooh_dir}", flush=True)
    outcar = Path(cwd) / "OUTCAR.RDip"
    del data["pq_index"]
    data["ads1"] = "non"  # type: ignore
    data["ads2"] = "OOH"  # type: ignore
    if os.path.exists("OUTCAR.RDip"):
        print(f"Adsorbate relaxation already exists in {ooh_dir}", flush=True)
        if not check_ase_database(
            "e_adsorbate_without_corrections", data, master=False
        ):
            print(f"Writing to local database for {ooh_dir}", flush=True)
            read_and_write_database(outcar, "e_adsorbate_without_corrections", data)
        control = True
    elif check_ase_database("e_adsorbate_without_corrections", data, master=True):
        print("In master database already", flush=True)
        control = True
    else:
        place_adsorbate(ooh_dir, data)
        control = relax_adsorbate(ooh_dir, data)
    os.chdir(cwd)
    if control:
        return True, data
    else:
        return False, None


if __name__ == "__main__":
    main()
