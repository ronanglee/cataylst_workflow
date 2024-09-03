from ase.io import read, write
from ase.db import connect
import os
import json


# arg1 = input("Input database name")
# arg2 = input("Input perqueue ID")

arg1 = 'pristine_vac'
db = connect(f'/home/energy/rogle/asm_orr_rxn/master_databases/{arg1}_master.db')

data = {'run_structure': '/home/energy/rogle/asm_orr_rxn/new_metals/cataylst_workflow/runs/structures/CoN4CZ+2N', 'name': 'CoN4CZ+2N', 'pristine': '/home/energy/rogle/asm_orr_rxn/new_metals/cataylst_workflow/runs/operating_stability/non_and_adsorbate/BEEF-vdW/CoN4CZ+2N/Co/1/non/non/non/non/non/non/non/non', 'metal': 'Co', 'pq_index': [2], 'dopant': '', 'adsorbate': '/home/energy/rogle/asm_orr_rxn/new_metals/cataylst_workflow/runs/operating_stability/non_and_adsorbate/BEEF-vdW/CoN4CZ+2N/Co/1/non/non/non/non/non/OOH/ontop_metal/1', 'carbon_structure': 'zigzag', 'base_dir': '/home/energy/rogle/asm_orr_rxn/new_metals/cataylst_workflow'}

file_name = '/home/energy/rogle/asm_orr_rxn/new_metals/cataylst_workflow/runs/operating_stability/non_and_adsorbate/BEEF-vdW/CoN4CZ+2N/Co/1/non/non/non/non/non/non/non/non/vac/vasp_rx'

del data['pq_index']
file = read(file_name + '/OUTCAR.RDip')

db.write(file, key_value_pairs={**data})
for row in db.select():
    print(row.name)

