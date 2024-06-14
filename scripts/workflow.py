import glob
import os
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup, DynamicWidthGroup # type: ignore

base_dir = Path(__file__).parent.parent

# Perqueue needs files not imported functions
syn_stab_run_e_xc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_xc.py"
syn_stab_run_e_c = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_c.py"
syn_stab_run_e_mxc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_mxc.py"
syn_stab_run_final = Path(base_dir) / "scripts" / "synthesis_stability" / "run_final.py"

testing = Path(base_dir) / "scripts" / "synthesis_stability" / "test.py"
testing1 = Path(base_dir) / "scripts" / "synthesis_stability" / "testing1.py"


syn_stab_run_e_m_on_c = (
    Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_m_on_c.py"
)

# TODO
# heat of formation
# Need to check to see percentage of each type of structure is left at the end of the workflow in comparison to the beginning # noqa
# Setting inital magnetic moments

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

    total_structures = dopant_structures*len(dopants) + template_structures
    return run_structures, total_structures

metals = ['Pt', 'Pd']
dopants = ["B", "S"] # TODO DOESNT WORK FOR MORE THAN 1 DOPANT!!!!!!!!!!!!!!!v
run_structures, num_total_structures = generate_input_files(metals, base_dir, dopants)

test_run_structs = []
for i in run_structures:
    if (
        # "PtN4C" == i.split("/")[-1]
        # or "PtN4CA" == i.split("/")[-1]
        # or 
        "PtN4CZ" == i.split("/")[-1]
        or "PtN3B1CZ+2N1B" == i.split("/")[-1]
    ):
        test_run_structs.append(i)

run_structures = test_run_structs

# Temporary for testing
count = 0
for i in run_structures:
    if 'B' in run_structures:
        count += len(dopants)
    else:
        count += 1
# run_structures = run_structures[0:6]

# resources = "56:1:xeon56:50h"
resources = "40:1:xeon40el8_768:10m"
# resources = '1:local:10m'

t1_data = {}
# for idx, e_c in enumerate(["zigzag", "armchair", "bulk"]):
for idx, e_c in enumerate(["zigzag"]):
    t1_data[idx] = {"base_dir": str(base_dir), "carbon_structure": e_c}


t2_data = {}
for idx1, e_c in enumerate(["zigzag"]):
    for idx2, metal in enumerate(metals):
            t2_data[idx1, idx2] = {"base_dir": str(base_dir), "carbon_structure": e_c, "metal": metal, 'all_run_structures': run_structures, 'dopants': dopants}
 
# for idx, data in enumerate(t1_data):
for idx, _ in enumerate(t1_data):
    t1 = t1_data[idx]
    t1 = Task(str(syn_stab_run_e_c), t1, resources=resources)
    extracted_data = [value for key, value in t2_data.items() if key[0] == idx]
    for data in extracted_data: # loop for each metal
        len_data = count
        t2 = Task(str(syn_stab_run_e_m_on_c), data, resources=resources) 
        # t22 = Task(str(check), data, resources=resources)
        t3 = Task(str(syn_stab_run_e_xc), {}, resources=resources)
        t4 = Task(str(syn_stab_run_e_mxc), {}, resources=resources)
        t5 = Task(str(syn_stab_run_final), {}, resources=resources)
        dwg1 = StaticWidthGroup(t3, width=len_data)
        dwg2 = StaticWidthGroup(t4, width=len_data)
        dwg3 = StaticWidthGroup(t5, width=len_data)
        swf1 = Workflow({t1: [], t2: [t1], dwg1: [t2], dwg2: [t2], dwg3: [t2, dwg1, dwg2]})
        with PersistentQueue() as pq:
            pq.submit(swf1)