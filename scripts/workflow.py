import glob
import os
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup, DynamicWidthGroup # type: ignore
from perqueue.constants import DYNAMICWIDTHGROUP_KEY

base_dir = Path(__file__).parent.parent

# Perqueue needs files not imported functions
syn_stab_run_e_xc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_xc.py"
syn_stab_run_e_c = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_c.py"
syn_stab_run_e_mxc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_mxc.py"
syn_stab_run_final = Path(base_dir) / "scripts" / "synthesis_stability" / "run_final.py"

testing = Path(base_dir) / "scripts" / "synthesis_stability" / "run_testing.py"
# testing1 = Path(base_dir) / "scripts" / "synthesis_stability" / "testing1.py"


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

    total_num_structures = dopant_structures*len(dopants) + template_structures
    return run_structures, total_num_structures

def gather_structs(data: dict) -> dict:
    '''Gather all the structures for next part of workflow.
    
    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        
    Returns:
        dict: Dictionary containing the run structure, base directory and metal.
    '''
    workflow_data = {}
    run_structures = data["all_run_structures"]
    metal = data["metal"]
    dopants = data["dopants"]
    carbon_structure = data["carbon_structure"]
    wanted_structures = []
    for struc_path in run_structures:
        structure = read(Path(struc_path) / Path("init.POSCAR"))
        struc_metal: str = 'Pt'
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
        if dopants and 'B' in structure.get_chemical_symbols():
            for dopant in dopants:
                data = {'base_dir': data['base_dir'], 'run_structure': str(struc_path.replace('B', dopant)), 'carbon_structure': carbon_structure, 'metal': metal, 'dopant': dopant}
                workflow_data[idx] = data
                idx += 1
        else:
            workflow_data[idx] = {'base_dir': data['base_dir'], 'run_structure': str(struc_path), 'carbon_structure': carbon_structure, 'metal': metal, 'dopant': ''}
            idx += 1

    return workflow_data


metals = ['Pt', 'Pd']
dopants = ["B", "S"] # TODO DOESNT WORK FOR MORE THAN 1 DOPANT!!!!!!!!!!!!!!!v
run_structures, num_total_structures = generate_input_files(metals, base_dir, dopants)

# Temporary for testing
count = 0
for i in run_structures:
    if 'B' in run_structures:
        count += len(dopants)
    else:
        count += 1
        
def gather_structs(data: dict) -> dict:
    '''Gather all the structures for next part of workflow.
    
    Args:
        data (dict): Dictionary containing the run structure, base directory and metal.
        
    Returns:
        dict: Dictionary containing the run structure, base directory and metal.
    '''
    workflow_data = {}
    run_structures = data["all_run_structures"]
    metal = data["metal"]
    dopants = data["dopants"]
    carbon_structure = data["carbon_structure"]
    wanted_structures = []
    for struc_path in run_structures:
        structure = read(Path(struc_path) / Path("init.POSCAR"))
        struc_metal: str = 'Pt'
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
        if dopants and 'B' in structure.get_chemical_symbols():
            for dopant in dopants:
                data = {'base_dir': data['base_dir'], 'run_structure': str(struc_path.replace('B', dopant).replace('Pt', metal)), 'carbon_structure': carbon_structure, 'metal': metal, 'dopant': dopant}
                workflow_data[idx] = data
                idx += 1
        else:
            workflow_data[idx] = {'base_dir': data['base_dir'], 'run_structure': str(struc_path).replace('Pt', metal), 'carbon_structure': carbon_structure, 'metal': metal, 'dopant': ''}
            idx += 1
    return workflow_data

# run_structures = run_structures[0:6]

# resources = "56:1:xeon56:50h"
# resources = "40:1:xeon40el8_768:10m"
resources = '1:local:10m'


test_run_structs = []
for i in run_structures:
    if (
        # "PtN4C" == i.split("/")[-1]
        "PtN4CA" == i.split("/")[-1]
        or 
        "PtN4CZ" == i.split("/")[-1]
        or "PtN3B1CZ+2N1B" == i.split("/")[-1]
    ):
        test_run_structs.append(i)

run_structures = test_run_structs


t1_data = {}
# for idx, e_c in enumerate(["zigzag", "armchair", "bulk"]):
for idx, e_c in enumerate(["zigzag", 'armchair']):
    t1_data[idx] = {"base_dir": str(base_dir), "carbon_structure": e_c}


t2_data = {}
idx = 0
for e_c in ["zigzag", 'armchair']:
    for metal in metals:
        t2_data[idx] = {"base_dir": str(base_dir), "carbon_structure": e_c, "metal": metal, 'all_run_structures': run_structures, 'dopants': dopants}
        idx +=1

all_structures = []
for key, value in t2_data.items():
    # print(gather_structs(value))
    # print('\n')
    all_structures.append(gather_structs(value))

t2s = []
for idx, t1_data in enumerate(t1_data.values()):
    t1 = Task(str(syn_stab_run_e_c), t1_data, resources=resources)
    # for strucs in all_structures[idx*len(metals)+1].values():
    for strucs in all_structures[idx*len(metals)].values():
        t2 = Task(str(syn_stab_run_e_xc), strucs, resources=resources)
        t2s.append(t2)
    t2 = Task(str(syn_stab_run_e_xc), all_structures[idx*len(metals)+1], resources=resources)
    #     t2s.append(t2)
    swg1 = StaticWidthGroup(t2, width=len(all_structures[idx*len(metals)]))
    # print(swg1)
    t2_wfs = []
    if idx == 0:
        f_idx = 0
    else:
        f_idx =1
    for idx2, t2_dat in enumerate(list(t2_data.values())[idx*(len(metals)):((idx)*len(metals) + len(metals))]):
        t3 = Task(str(syn_stab_run_e_m_on_c), t2_dat, resources=resources)
        t4 = Task(str(syn_stab_run_e_mxc), {}, resources=resources)
        t5 = Task(str(syn_stab_run_final), {}, resources=resources)
        # wff = Workflow({items: [t1] for items in t2s})
        swg3 = StaticWidthGroup([t4, t5], width=len(all_structures[idx*len(metals)]))
        swf = Workflow({t3: [t1], swg3: [t3, swg1]})
                
        t2_wfs.append(swf)
        # print(swf)
        # exit()
    wff = Workflow({items: [] for items in t2_wfs})
    wf = Workflow({t1: [], swg1: [t1], wff: []})
    with PersistentQueue() as pq:
        pq.submit(wf)
        




    # for all_struct in all_structures:
#     t5 = Task(str(syn_stab_run_e_mxc), all_structures, resources=resources)
#     t6 = Task(str(syn_stab_run_final), all_structures, resources=resources)
#     swg2 = StaticWidthGroup([t5, t6], width=len(all_structures))

# wf = Workflow({t2: [], swg1: [t2], t3: [t2], swg2: [t4]})
# with PersistentQueue() as pq:
#     pq.submit(wf)
# t2 = Task(str(syn_stab_run_e_m_on_c), t22_data, resources=resources)
# swg2 = StaticWidthGroup(t2, width=len(t22_data))
# print(swg2)

# for value1, metal in zip(t22_data.values(), metals):
#     print(value1)
#     t2 = Task(str(syn_stab_run_e_m_on_c), value1, resources=resources)
#     t3_items = []
#     t4_items = []
#     t5_items = []
#     valllue1 = value1
#     valllue1['all_run_structures'] = run_structures
#     valllue1['dopants'] = dopants
#     workflow_data = gather_structs(valllue1)
#     # print(workflow_data)
#     # exit()
#     new_data = {}
#     idx = 0
#     for value2 in workflow_data.values():
#         if value2['metal'] == metal:
#             new_data[idx] = value2
#             idx += 1
#     print(new_data)
#     t3 = Task(str(syn_stab_run_e_xc), {}, resources=resources)
#     t4 = Task(str(syn_stab_run_e_mxc), {}, resources=resources)
#     t5 = Task(str(syn_stab_run_final), {}, resources=resources)
#     t3_items.append(t3)
#     swg = DynamicWidthGroup([t3, t4])
#     wf = Workflow({t2: [], swg: [t2]})
#     # wf = Workflow({t1: []})
#     with PersistentQueue() as pq:
#         pq.submit(wf)


# t1 = Task(str(syn_stab_run_e_c), t1_data, resources=resources)

# t2 = Task(str(syn_stab_run_e_m_on_c), t22_data, resources=resources)
# swg = StaticWidthGroup(t2, width=len(t22_data))
# t3 = Task(str(syn_stab_run_e_xc), {}, resources=resources)
# swg2 = StaticWidthGroup(t3, width=len(run_structures))


# swg3 = StaticWidthGroup(t2, width=len(structures[0]))


# t4 = Task(str(syn_stab_run_e_mxc), {}, resources=resources)
# swg4 = StaticWidthGroup(t4, width=len(all_structures))
# t5 = Task(str(syn_stab_run_final), {}, resources=resources)
# swg5 = StaticWidthGroup(t5, width=len(all_structures))

# # wf = Workflow({swg1: [], swg2: [swg1]})
# wf = Workflow({swg1: [], swg2: [swg1], swg3: [swg1], swg4: [swg2], swg5: [swg2]})

# # wf = Workflow({swg3: []})
# # t1 = Task(str(syn_stab_run_e_c), t1_data, resources=resources)
# # swg1 = StaticWidthGroup(t1, width=len(t1_data))

