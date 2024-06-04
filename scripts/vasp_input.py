import os


def vasp_input():
    params = {
        "encut": 600,  # plane-wave cutoff in eV
        "kpts": (3, 3, 1),  # k-point sampling
        # smearing method: -1 for Fermi, 0 for Gaussian, see documentation for other
        "ismear": -1,
        "sigma": 0.1,  # width of the smearing in eV
        "ispin": 2,  # spin polarized? (1: no, 2: yes)
        "lorbit": 11,  # print magmoms to the OUTCAR
        # use 'Fast' in case of convergence problems try 'Normal' or 'All' (see documentation)
        "algo": "Normal",
        "ediff": 1e-05,  # break condition for the electronic SC-loop
        "nelmin": 0,  # Force the minimum number of SC-steps
        # break condition for ionic relaxation, positive: energy; negative: forces
        "ediffg": -0.02,
        # relaxation scheme (0: relax ions only; see documentation for other)
        "isif": 4,
        # Use of symmetry. 0: switch off symmetry. Default: 2 for PAW (1 for US-PP)
        "isym": 2,
        "prec": "Accurate",  # use 'Accurate' (see documentation)
        "nelm": 3000,  # Max. number of SC steps (default: 60)
        # Scaling constant for the forces (or time step for ibrion=0)
        "potim": 0.1,
        "lcharg": True,  # write charge densities to CHGCAR file?
        # construct 'initial' chg. dens.: 0: from wavefun., 1: from CHGCAR, 2: superpos'n of atomic
        "icharg": 2,
        "lwave": False,  # write wavefunctions it WAVECAR file?
        "lcorr": True,  # Corrections to the forces
        "lmaxmix": 4,  # See manual (or rather VASP wiki)
        "lasph": True,  # See manual
        "lreal": "Auto",  # See manual
        # See manual. Shouldn't affect energy, only performance. 4 is default.
        "nsim": 4,
        # See manual. True is good for metals and small gap semiconductors.
        "ldiag": True,
        "xc": "PBE",
        # 'gga'     : 'bf',
        # 'zab_vdw' : -1.8867,
        # 'luse_vdw': True,
        # 'maxmix'  : 200,
        # 'ncore'   : 40
    }
    return params


def synthesis_stability_run_vasp(directory: os.PathLike, vasp_parameters: dict) -> None:
    """Run the VASP calculations.

    Args:
        directory (Path): Path to the directory where to find the initial geometry.
        vasp_parameters (dict): Dictionary containing the VASP parameters.

    Returns:
        None
    """
    os.chdir(directory)
    # Heine's standard mixing parameters
    vasp_parameters["amin"] = 0.1
    vasp_parameters["amix"] = 0.02
    vasp_parameters["bmix"] = 1.0
    if vasp_parameters["ispin"] == 2:  # Check if spin polarized
        vasp_parameters["amix_mag"] = 0.08
        vasp_parameters["bmix_mag"] = 1.0

    atoms = read(directory / "init.POSCAR")
    # Begin with static calculation without dipole correction
    vasp_parameters["nelmdl"] = -8
    vasp_parameters["ibrion"] = 2
    vasp_parameters["nsw"] = 0
    paramscopy = vasp_parameters.copy()
    calc = Vasp(**paramscopy)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
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
    vasp_parameters["nsw"] = 999
    vasp_parameters["lwave"] = True
    paramscopy = params.copy()
    calc = Vasp(**paramscopy)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
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
