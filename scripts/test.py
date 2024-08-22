from pathlib import Path

data = {'base_dir': '/home/energy/rogle/asm_orr_rxn/catalyst_workflow', 'run_structure': '/home/energy/rogle/asm_orr_rxn/catalyst_workflow/runs/structures/CoN4CZ+1B1S', 'carbon_structure': 'zigzag', 'dopant': 'SB', 'run_dir': '/home/energy/rogle/asm_orr_rxn/catalyst_workflow/runs/operating_stability/e_xch_implicit_solv/MN4CZ+1B1S/1H_config1'}

database: dict = {}
struc_name = Path(data["run_structure"]).stem
config = Path(data['run_dir']).stem
database[struc_name] = {config: {'xch_corrections': 'ronan'}}


corrected = database[struc_name][config]['xch_corrections']
print(corrected)
