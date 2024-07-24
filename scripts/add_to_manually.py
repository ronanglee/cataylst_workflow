from ase.db import connect
from ase.io import read

data =  {'run_structure': '/home/energy/rogle/asm_orr_rxn/catalyst_workflow/runs/structures/CoN4C+1B_second_coordination', 'ads1': 'non', 'base_dir': '/home/energy/rogle/asm_orr_rxn/catalyst_workflow', 'dopant': 'B', 'name': 'CoN4C+1B_second_coordination', 'metal': 'Co', 'carbon_structure': 'bulk', 'pq_index': [1], 'ads2': 'OOH'}
db = connect('/home/energy/rogle/asm_orr_rxn/catalyst_workflow/runs/databases/pristine_implicit.db')
outcar = read('OUTCAR')
del data['pq_index']
db.write(outcar, key_value_pairs={**data})

