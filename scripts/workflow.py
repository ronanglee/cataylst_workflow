import glob
import os
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import StaticWidthGroup, Workflow  # type: ignore
from utils import gather_structs, run_logger  # type: ignore

base_dir = Path(__file__).parent.parent

# Perqueue needs files not imported functions
syn_stab_run_e_xc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_xc.py"
syn_stab_run_e_c = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_c.py"
syn_stab_run_e_mxc = Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_mxc.py"
syn_stab_run_final = Path(base_dir) / "scripts" / "synthesis_stability" / "run_final.py"
syn_stab_run_e_m_on_c = (
    Path(base_dir) / "scripts" / "synthesis_stability" / "run_e_m_on_c.py"
)

oper_stab_run_ooh = Path(base_dir) / "scripts" / "operating_run_scripts" / "run_ooh.py"
oper_stab_run_pristine = (
    Path(base_dir) / "scripts" / "operating_run_scripts" / "run_pristine.py"
)
oper_stab_run_ads = (
    Path(base_dir) / "scripts" / "operating_run_scripts" / "run_adsorbate.py"
)
oper_stab_run_vib = (
    Path(base_dir) / "scripts" / "operating_run_scripts" / "run_vibration.py"
)
oper_stab_run_xch = Path(base_dir) / "scripts" / "operating_run_scripts" / "run_xch.py"

# TODO
# heat of formation
# Need to check to see percentage of each type of structure is left at the end of the workflow in comparison to the beginning # noqa
# Setting inital magnetic moments???????
# Need to implement a logger for each step where we filter out the structures that are not stable # noqa
# Raul variational autoencoder

if os.path.exists("workflow.log"):
    os.remove("workflow.log")

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


metals = ["Pt", "Pd"]
dopants = ["B", "S"]
run_structures = generate_input_files(metals, base_dir, dopants)

# resources = "56:1:xeon56:50h"
# resources = "40:1:xeon40el8_768:10m"
# resources = '1:local:10m'
resources = "5:1:xeon24el8_test:10m"
# resources =  '40:1:xeon40el8:10m'

run_logger("Starting workflow", str(__file__), "info")
# Should have an end workflow somewhere
# exit()

test_run_structs = []
for i in run_structures:
    if (
        # "PtN4C" == i.split("/")[-1]
        "PtN4CA" == i.split("/")[-1]
        or "PtN4CZ" == i.split("/")[-1]
        or "PtN3B1CZ+2N1B" == i.split("/")[-1]
    ):
        test_run_structs.append(i)
run_structures = test_run_structs

# Warnings
print("Dont have good vasp parameters")
print("ARE THESE STRUCTURES EVEN OPTIMISED, especially the adsorbate structures???")
print("no false return on functions")
print("no ionic/scf checks")
print("end of pristine run incorrect")
print("end of adsorption run incorrect")
print("Need to check hydrogens most stable structure")

t1_data = {}
# for idx, e_c in enumerate(["zigzag", "armchair", "bulk"]):
for idx, e_c in enumerate(["zigzag", "armchair"]):
    t1_data[idx] = {"base_dir": str(base_dir), "carbon_structure": e_c}

t2_data = {}
idx = 0
for e_c in ["zigzag", "armchair"]:
    for metal in metals:
        t2_data[idx] = {
            "base_dir": str(base_dir),
            "carbon_structure": e_c,
            "metal": metal,
            "all_run_structures": run_structures,
            "dopants": dopants,
        }
        idx += 1

all_structures = []
for key, value in t2_data.items():
    all_structures.append(gather_structs(value))

for idx, t1_data in enumerate(t1_data.values()):  # type: ignore
    t1 = Task(str(syn_stab_run_e_c), t1_data, resources=resources)
    t2 = Task(
        str(syn_stab_run_e_xc),
        all_structures[idx * len(metals) + 1],
        resources=resources,
    )
    t2h = Task(
        str(oper_stab_run_xch),
        all_structures[idx * len(metals) + 1],
        resources=resources,
    )

    swg1 = StaticWidthGroup([t2, t2h], width=len(all_structures[idx * len(metals)]))
    # swg1 = StaticWidthGroup([t2], width=len(all_structures[idx*len(metals)]))
    if idx == 0:
        f_idx = 0
    else:
        f_idx = 1
    t2_wfs = []
    for t2_dat in list(t2_data.values())[
        idx * (len(metals)) : ((idx) * len(metals) + len(metals))  # noqa
    ]:
        t3 = Task(str(syn_stab_run_e_m_on_c), t2_dat, resources=resources)
        t4 = Task(str(syn_stab_run_e_mxc), {}, resources=resources)
        t5 = Task(str(syn_stab_run_final), {}, resources=resources)
        t6 = Task(str(oper_stab_run_ooh), {}, resources=resources)
        t7 = Task(str(oper_stab_run_pristine), {}, resources=resources)
        t8 = Task(str(oper_stab_run_ads), {}, resources=resources)
        t9 = Task(str(oper_stab_run_vib), {}, resources=resources)
        swg3 = StaticWidthGroup(
            [t4, t5, t6, t7, t8, t9], width=len(all_structures[idx * len(metals)])
        )
        swf = Workflow({t3: [t1], swg3: [t3]})
        t2_wfs.append(swf)
    swf2 = Workflow({items: [] for items in t2_wfs})
    wf = Workflow({swg1: [], t1: [], swf2: []})
    with PersistentQueue() as pq:
        pq.submit(wf)
