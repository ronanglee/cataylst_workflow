def vasp_input() -> dict:
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
        "ediff": 1e-06,  # break condition for the electronic SC-loop
        # break condition for ionic relaxation, positive: energy; negative: forces
        "ediffg": -0.025,
        "nelmin": 0,  # Force the minimum number of SC-steps
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
        "symprec": 1e-08,  # See manual
        # See manual. Shouldn't affect energy, only performance. 4 is default.
        "nsim": 4,
        # See manual. True is good for metals and small gap semiconductors.
        "ldiag": True,
        "xc": "PBE",
        "gga": "bf",
        "amin": 0.1,
        "amix": 0.02,
        "bmix": 1.0,
        "amix_mag": 0.08,
        "bmix_mag": 1.0,
        "ibrion": 2,
        "zab_vdw": -1.8867,
        "luse_vdw": True,
        "maxmix": 200,
        "ncore": 4,
    }
    return params
