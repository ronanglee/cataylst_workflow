import glob
import os
from pathlib import Path

from ase.io import read, write  # type: ignore
from perqueue import PersistentQueue  # type: ignore
from perqueue.task_classes.task import Task  # type: ignore
from perqueue.task_classes.task_groups import StaticWidthGroup, Workflow  # type: ignore
from utils import gather_structs, get_xch_structs, run_logger  # type: ignore

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
oper_stab_run_xch_solv_implicit = (
    Path(base_dir) / "scripts" / "operating_run_scripts" / "run_xch_solv_implicit.py"
)
oper_stab_run_xc_solv_implicit = (
    Path(base_dir) / "scripts" / "operating_run_scripts" / "run_xc_solv_implicit.py"
)
relative_stab = Path(base_dir) / "scripts" / "operating_stability" / "run_e_r.py"
act_sel = Path(base_dir) / "scripts" / "activity_selectivity" / "run_act_sel.py"


cwd = Path(__file__).parent
run_folder = Path(base_dir) / "runs"
data_base_folder = run_folder / "databases"
run_folder.mkdir(exist_ok=True)
data_base_folder.mkdir(exist_ok=True)
if not os.path.exists(base_dir / "template_structures"):
    raise FileNotFoundError(
        "Directory not found: template_structures. Need to have template_structures folder in the base directory."
    )


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
        run_folder = Path(base_dir) / "runs" / "structures"
        run_structures = []
        # LARGE SPACES ARE FOR THE PATCHY CODE
        if ["SB"] == dopants:
            dop_glob = "*S*B"
        elif ["O"] == dopants:
            dop_glob = "*O*"
        else:
            dop_glob = "*"
        dopant_files = glob.glob(
            str(Path(base_dir) / "template_structures" / "test_structs" / "dopant")
            + "/"
            + f"{dop_glob}"
            + "/POSCAR.opt"
        )
        dopant_structures = [f for f in dopant_files]
        for idx, input in enumerate(dopant_structures):
            structure = read(input)
            for metal in metals:
                copy_structure = structure.copy()
                for atom in copy_structure:
                    if atom.symbol == "Pt":
                        atom.symbol = metal
                    save_dir = run_folder / dopant_files[idx].split("/")[-2].replace(
                        "Pt", metal
                    )
                run_structures.append(str(save_dir))
                save_dir.mkdir(exist_ok=True)
                write(save_dir / "init.POSCAR", copy_structure)

                new_dopants = []
                for dopant in dopants:
                    if dopant != "O":
                        if dopant != "SB":
                            new_dopants.append(dopant)
                if "S" not in copy_structure.get_chemical_symbols():
                    if "O" not in copy_structure.get_chemical_symbols():
                        for dopant in new_dopants:
                            copy_structure = structure.copy()
                            for atom in copy_structure:
                                if atom.symbol == "Pt":
                                    atom.symbol = metal
                                if atom.symbol == "B":
                                    atom.symbol = dopant
                            save_dir = run_folder / dopant_files[idx].split("/")[
                                -2
                            ].replace("Pt", metal).replace("B", dopant)
                            # run_structures.append(str(save_dir))
                            save_dir.mkdir(exist_ok=True)
                            write(save_dir / "init.POSCAR", copy_structure)

    return run_structures


metals = ["Co" "Pt"]
dopants = ["S", "B", "O", "SB"]
run_structures = generate_input_files(metals, base_dir, dopants)

xeon40 = "40:1:xeon40el8:25h"

run_logger("Starting workflow", str(__file__), "info")

t1_data = {}
for idx, e_c in enumerate(["bulk", "zigzag", "armchair"]):
    t1_data[idx] = {"base_dir": str(base_dir), "carbon_structure": e_c}

t2_data = {}
idx = 0
for idx, e_c in enumerate(["bulk", "zigzag", "armchair"]):
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
    all_structures.append(gather_structs(value, metals))

zigzag_count = 0
bulk_count = 0
armchair_count = 0
run_structs_logger = []
for outer_dict in all_structures:
    for key, inner_dict in outer_dict.items():
        if inner_dict.get("carbon_structure") == "zigzag":
            zigzag_count += 1
        elif inner_dict.get("carbon_structure") == "bulk":
            bulk_count += 1
        elif inner_dict.get("carbon_structure") == "armchair":
            armchair_count += 1
        inner_struc = str(Path(inner_dict.get("run_structure")).stem)
        run_structs_logger.append(inner_struc)

total_count = zigzag_count + bulk_count + armchair_count
run_logger(f"Starting structures {run_structs_logger}", str(__file__), "info")
run_logger(
    f"Starting counts - bulk count {{{bulk_count}}}, zigzag count {{{zigzag_count}}}, armchair count {{{armchair_count}}} - Total {{{total_count}}}",  # noqa
    str(__file__),
    "info",
)
run_logger(f"metals {metals}, dopants {dopants}", str(__file__), "info")


for idx, t1_data in enumerate(t1_data.values()):  # type: ignore
    t1 = Task(str(syn_stab_run_e_c), t1_data, resources=xeon40)
    t2_wfs = []
    for structure in all_structures[idx * len(metals)].values():
        t2 = Task(
            str(syn_stab_run_e_xc),
            structure,
            resources=xeon40,
        )
        tsub2 = Task(
            str(oper_stab_run_xc_solv_implicit),
            structure,
            resources=xeon40,
        )

        e_xch_data = get_xch_structs(structure)
        t2h = Task(str(oper_stab_run_xch_solv_implicit), e_xch_data, resources=xeon40)

        swg1 = StaticWidthGroup([t2h], width=len(e_xch_data))
        wf = Workflow({tsub2: [t2], swg1: [tsub2]})
        t2_wfs.append(wf)

    t3_wfs = []
    for t2_dat in list(t2_data.values())[
        idx * (len(metals)) : ((idx) * len(metals) + len(metals))  # noqa
    ]:
        t2_dat["metals"] = metals
        t3 = Task(str(syn_stab_run_e_m_on_c), t2_dat, resources=xeon40)
        t4 = Task(str(syn_stab_run_e_mxc), {}, resources=xeon40)
        t5 = Task(str(syn_stab_run_final), {}, resources="1:1:xeon24el8:10m")
        t6 = Task(str(oper_stab_run_ooh), {}, resources="1:1:xeon24el8:10m")
        t7 = Task(str(oper_stab_run_pristine), {}, resources=xeon40)
        t8 = Task(str(oper_stab_run_ads), {}, resources=xeon40)
        t9 = Task(str(oper_stab_run_vib), {}, resources=xeon40)
        t10 = Task(str(relative_stab), {}, resources="1:1:xeon24el8:10m")
        t11 = Task(str(act_sel), {}, resources="1:1:xeon24el8:10m")
        swg2 = StaticWidthGroup(
            [t4, t5, t6, t7, t8, t9, t10, t11],
            width=len(all_structures[idx * len(metals)]),
        )
        swf = Workflow({t3: [t1], swg2: [t3]})
        t3_wfs.append(swf)

    swf2 = Workflow({items: [] for items in t2_wfs})
    swf3 = Workflow({items: [] for items in t3_wfs})
    wf = Workflow({t1: [], swf2: [], swf3: []})
    with PersistentQueue() as pq:
        pq.submit(wf)
