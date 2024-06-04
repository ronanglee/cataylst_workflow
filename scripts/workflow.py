from synthesis_stability.generate_run import main
import glob
import os
import sys
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
# type: ignore

base_dir = Path(__file__).parent.parent
cwd = Path(__file__).parent
run_folder = Path(base_dir) / "runs"
run_folder.mkdir(exist_ok=True)


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
    # print(f"Calculating for {len(run_files) *
    #       len(metals) * len(dopant)} structures")
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

synthesis_stability = (
    Path(base_dir) / "scripts" / "synthesis_stability" / "generate_run.py"
)

sys.path.append("/home/energy/rogle/asm_orr_rxn/catalyst_workflow/scripts")

for struc_path in run_structures:
    structure = read(Path(struc_path) / Path("init.POSCAR"))
    dopant = None
    struc_metal = None
    for symbol in structure.get_chemical_symbols():
        if symbol in metals:
            struc_metal = symbol
        if dopants:
            if symbol in dopants:
                dopant = symbol
    if structure[35].symbol == struc_metal:
        carbon_structure = "armchair"
    elif structure[45].symbol == struc_metal or structure[38].symbol == struc_metal:
        carbon_structure = "zigzag"
    elif structure[62].symbol == struc_metal:
        carbon_structure = "bulk"
    data = {
        "base_dir": str(base_dir),
        "run_structure": struc_path,
        "carbon_structure": carbon_structure,
        "metal": struc_metal,
        "dopant": dopant,
    }
    # print(data)
    # main(**data)
    # print("RUNNING WITHOUT VDW CORRECTION")

    t1 = Task(str(synthesis_stability), data, "24:1:xeon24:1m")
    wf = Workflow({t1: []})
    with PersistentQueue() as pq:
        pq.submit(wf)

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
