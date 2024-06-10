import glob
import os
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import Workflow  # type: ignore

base_dir = Path(__file__).parent.parent

# Perqueue needs files not imported functions
syn_stab_run_e_xc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_xc.py"
syn_stab_run_e_c = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_c.py"
syn_stab_run_e_mnx = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_mxc.py"
syn_stab_run_e_m_on_c = (
    Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_m_on_c.py"
)

# TODO
# heat of formation
# Need to check to see percentage of each type of structure is left at the end of the workflow in comparison to the beginning # noqa


try:  # type: ignore
    os.system(f"rm -r {base_dir}/runs")  # type: ignore
except:  # type: ignore # noqa
    pass  # type: ignore # noqa
cwd = Path(__file__).parent
run_folder = Path(base_dir) / "runs"
data_base_folder = run_folder / "databases"
run_folder.mkdir(exist_ok=True)
data_base_folder.mkdir(exist_ok=True)


def generate_input_files(
    metals: list, base_dir: os.PathLike, dopants: list | None = None
) -> list:
    """Generates input files for the calculations. Default metal is Pt and default dopant is B.
    All stuctures contain N as default. Template files found in template_structures folder.
    Need to have Pt and N (and B if want to add dopant) in the template file name.

    Args:
        metals (List): List of metals.
        dopants (List): List of dopants.

    Returns:
        run_structures (List): List of paths to the generated input files.
    """
    template_files = glob.glob(
        str(Path(base_dir) / "template_structures" / "test_structs") + "/Pt*/POSCAR.opt"
    )
    template_structures = [read(f) for f in template_files]
    run_folder = Path(base_dir) / "runs" / "structures"
    run_folder.mkdir(exist_ok=True)
    run_structures = []
    for idx, structure in enumerate(template_structures):
        for metal in metals:
            copy_structure = structure.copy()
            for atom in copy_structure:
                if atom.symbol == "Pt":
                    atom.symbol = metal
            save_dir = run_folder / template_files[idx].split("/")[-2].replace(
                "Pt", metal
            )
            run_structures.append(str(save_dir))
            save_dir.mkdir(exist_ok=True)
            write(save_dir / "init.POSCAR", copy_structure)
    if dopants:
        dopant_files = glob.glob(
            str(Path(base_dir) / "template_structures" / "test_structs" / "dopant")
            + "/*B*/POSCAR.opt"
        )
        dopant_structures = [read(f) for f in dopant_files]
        for idx, structure in enumerate(dopant_structures):
            for metal in metals:
                for dope in dopants:
                    copy_structure = structure.copy()
                    for atom in copy_structure:
                        if atom.symbol == "Pt":
                            atom.symbol = metal
                        if atom.symbol == "B":
                            atom.symbol = dope
                    save_dir = run_folder / dopant_files[idx].split("/")[-2].replace(
                        "Pt", metal
                    ).replace("B", dope)
                    run_structures.append(str(save_dir))
                    save_dir.mkdir(exist_ok=True)
                    write(save_dir / "init.POSCAR", copy_structure)

    return run_structures


metals = ["Pd", "Pt"]
dopants = ["B", "S"]
run_structures = generate_input_files(metals, base_dir, dopants)

test_run_structs = []
for i in run_structures:
    if (
        "PtN4C" == i.split("/")[-1]
        or "PtN4CA" == i.split("/")[-1]
        or "PtN4CZ" == i.split("/")[-1]
        or "PtN3B1CZ+2N1B" == i.split("/")[-1]
    ):
        test_run_structs.append(i)

run_structures = test_run_structs
# run_structures = run_structures[0:6]

resources = "24:1:xeon24:20h"

for e_c in ["armchair", "zigzag", "bulk"]:
    t1 = Task(
        str(syn_stab_run_e_c),
        {"base_dir": str(base_dir), "carbon_structure": e_c},
        resources,
    )
    t2_structs = []
    for metal in metals:
        t2 = Task(
            str(syn_stab_run_e_m_on_c),
            {"base_dir": str(base_dir), "carbon_structure": e_c, "metal": metal},
            resources,
        )
        t2_structs.append(t2)

    t3_structs = []
    t4_structs = []
    for struc_path in run_structures:
        structure = read(Path(struc_path) / Path("init.POSCAR"))
        struc_metal: str | None = None
        carbon_structure: str | None = None
        if metal in structure.get_chemical_symbols():
            struc_metal = metal
        else:
            continue
        if structure[35].symbol == struc_metal:
            carbon_structure = "armchair"
        elif structure[45].symbol == struc_metal or structure[38].symbol == struc_metal:
            carbon_structure = "zigzag"
        elif structure[62].symbol == struc_metal:
            carbon_structure = "bulk"
        if struc_metal == metal and carbon_structure == e_c:
            t3 = Task(
                str(syn_stab_run_e_xc),
                {"base_dir": str(base_dir), "run_structure": Path(struc_path)},
                resources,
            )
            t4 = Task(
                str(syn_stab_run_e_mnx),
                {"base_dir": str(base_dir), "run_structure": Path(struc_path)},
                resources,
            )
            t3_structs.append(t3)
            t4_structs.append(t4)
    t2s = {items: [t1] for items in t2_structs}
    swf1 = Workflow(t2s)
    t3s = {items: [swf1] for items in t3_structs}
    t4s = {items: [swf1] for items in t4_structs}
    wf1 = Workflow(t3s)
    wf2 = Workflow(t4s)
    wf = Workflow({wf1: [], wf2: []})
    with PersistentQueue() as pq:
        pq.submit(wf)


exit()
# for struc_path in run_structures:
#     structure = read(Path(struc_path) / Path("init.POSCAR"))
#     dopant: str | None = None
#     struc_metal: str | None = None
#     carbon_structure: str | None = None
#     for symbol in structure.get_chemical_symbols():
#         if symbol in metals:
#             struc_metal = symbol
#         if dopants:
#             if symbol in dopants:
#                 dopant = symbol
#     if structure[35].symbol == struc_metal:
#         carbon_structure = "armchair"
#     elif structure[45].symbol == struc_metal or structure[38].symbol == struc_metal:
#         carbon_structure = "zigzag"
#     elif structure[62].symbol == struc_metal:
#         carbon_structure = "bulk"
#     if not struc_metal:
#         print("No metal found in structure")
#     if not carbon_structure:
#         print("No carbon structure found in structure")
#     if not dopant:
#         print("No dopant found in structure")
#     data = {
#         "base_dir": str(base_dir),
#         "run_structure": struc_path,
#         "carbon_structure": carbon_structure,
#         "metal": struc_metal,
#         "dopant": dopant,
#     }
#     # main(**data)
#     t4 = Task(str(syn_stab_run_e_c), data, resources)
#     wf = Workflow({t4: []})
#     with PersistentQueue() as pq:
#         pq.submit(wf)
#
#     # print(data)
#     # main(**data)
#     # print("RUNNING WITHOUT VDW CORRECTION")
#
#     t3 = Task(str(syn_stab_run_e_xc), data, resources)
#     t4 = Task(str(syn_stab_run_e_mnx), data, resources)
# e_c_dir = base_dir / "runs" / "synthesis_stability" / "e_c" / f"{carbon_structure}"
# e_m_on_c_dir = (
#     base_dir
#     / "runs"
#     / "synthesis_stability"
#     / "e_m_on_c"
#     / f"{carbon_structure}_0N_0H"
#     / f"{struc_metal}"
# )
# #
# #     e_c_dir.exists() and e_m_on_c_dir.exists()
# # ):  # Skip as carbon structure already exists
# #     wf = Workflow({t1: [], t2: []})
# # elif e_m_on_c_dir.exists() and not e_c_dir.exists():
# #     e_c_dir.mkdir(exist_ok=True, parents=True)
# #     t3 = Task(str(syn_stab_run_e_c), data, resources)
# #     wf = Workflow({t1: [], t2: [], t3: []})
# # elif not e_m_on_c_dir.exists() and e_c_dir.exists():
# #     e_m_on_c_dir.mkdir(exist_ok=True, parents=True)
# #     t4 = Task(str(syn_stab_run_e_m_on_c), data, resources)
# #     wf = Workflow({t1: [], t2: [], t4: []})
# # elif not e_m_on_c_dir.exists() and not e_c_dir.exists():
# #     e_c_dir.mkdir(exist_ok=True, parents=True)
# #     e_m_on_c_dir.mkdir(exist_ok=True, parents=True)
# t3 = Task(str(syn_stab_run_e_c), data, resources)
# t4 = Task(str(syn_stab_run_e_m_on_c), data, resources)
# wf = Workflow({t1: [t3, t4], t2: [t3, t4]})
# # I should run the dependencies first and then the main task
# with PersistentQueue() as pq:
#     pq.submit(wf)

#
# for struc in run_structures:
#     data = {
#         "base_dir": str(base_dir / "runs"),
#         "run_structure": struc,
#         "metals": metals,
#         "dopant": dopant,
#     }
#     # '48:xeon24el8_test:10m')
#     t1 = Task(str(synthesis_stability), data, "24:1:xeon24el8_test:10s")
#     # t1_swg = StaticWidthGroup(t1, width=len(run_structures))
#     # # The first time, it runs for all the displacements.
#     # # It will generate multiple processes due to the parallelization.
#     # # I could run them in serialized mode, but it would take a long time.
#     # # Even when it marks an error at the end, the process is still running.
#     # t2 = Task(cwd/'supercell.py', args, '1:xeon24el8:50h') # Parallelization must be disabled...
#     # # Use 1 core xeon24 if not, then ask for a full node, if not, then use high-memory nodes.
#     # t3 = Task(cwd/'phonons.py', args, '1:xeon24el8:1h')
#     # t4 = Task(cwd/'scf.py',args, '1:xeon24el8:50h') #'40:xeon40el8:1h') # Done quickly
#     # t5 = Task(cwd/'dipolemoment.py', args, '1:xeon24el8:10h')
#     # t6 = Task(cwd/'gmatrix.py', args, '1:xeon24el8:1h') # Not parallelizable
#     # t7 = Task(cwd/'raman.py', args, '1:xeon24el8:1h') # Run locally.
#     # t8 = Task(cwd/'plot_spectrum.py', args, '1:xeon24el8:1h')
#
#     # swf1 = Workflow({t2: [], t3: []})
#     # swf2 = Workflow({t4: [], t5: [t4]})
#     # wf = Workflow({t1: [], swf1: [t1], t6: [swf1], t7: [t6, swf2], t8: [t7]})
#     wf = Workflow({t1: []})
#
#     with PersistentQueue() as pq:
#         pq.submit(wf)
# #     #pq.submit(t1)
