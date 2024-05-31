from ase.io import read, write
import os
import glob
from pathlib import Path
from collections import Counter

files = glob.glob('/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/*/POSCAR.opt')
tia_path = '/home/energy/tipapa/Yang_work/Pt/'
coordination_sphere_1 = ['3', '4'] # everything else has weird pores in them
coordination_sphere_2 = []

n_counter = []

def look_at_data():
    print(files)
    for file in files:
        # print(file)
        structure = read(file)
        m_idx = [i for i, x in enumerate(structure.get_chemical_symbols()) if x == "Pt"]  
        n_idx = [i for i, x in enumerate(structure.get_chemical_symbols()) if x == "N"]
        distances = []
        for i in m_idx:
            for j in n_idx:
                distances.append(structure.get_distance(i, j))

        sphere_one_ount = 0
        sphere_two_count = 0
        sphere_three_count = 0
        for i in sorted(distances):
            if i < 2.6:
                sphere_one_ount += 1
            elif i < 3.6 and i > 2.6:
                sphere_two_count += 1
            elif i > 3.6:
                sphere_three_count += 1
        co_str = f'N{sphere_one_ount}+{sphere_two_count}N+{sphere_three_count}N'
        n_count = sphere_one_ount + sphere_two_count + sphere_three_count
        n_counter.append(n_count)
        print(file, co_str, n_count)
            
    print(Counter(n_counter))
    
look_at_data()

def get_structures():
    for file in files:
        if 'png' in file:
            dir = file.split('/')[-1].replace('.png', '')
            try:
                os.mkdir(dir)
            except:
                pass
            structure_path = Path(tia_path) / file.split('/')[-1].replace('.png', '') / 'POSCAR.opt'
            os.system(f'cp {structure_path} {dir}')
            structure = read(structure_path)
            write(f'./{dir}/{file.split('/')[-1]}', structure)
 
# get_structures()
           
# def copy_structs():
#     files = glob.glob('/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/PtN4CZ*')
#     print(files)
    
# copy_structs()