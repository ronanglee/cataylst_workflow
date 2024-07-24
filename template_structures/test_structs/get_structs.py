import glob
import os
from collections import Counter
from pathlib import Path
from ase.visualize import view
import re
from ase.io import read, write

files = glob.glob(
    "/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/*/POSCAR.opt"
)

parent_dir = "/home/energy/rogle/asm_orr_rxn/catalyst_workflow/template_structures/test_structs/"
files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
print(len(files))
# files = [file for file in files if "B" not in file]
# files = files[32:] + ['/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/PtN3CZ/POSCAR.opt', '/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/PtN4C/POSCAR.opt', '/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/PtN4CA/POSCAR.opt', '/home/energy/rogle/asm_orr_rxn/workflow/template_structures/test_structs/PtN4CZ/POSCAR.opt']
# for idx, file in enumerate(files):
#     print(idx, file)
# exit()



tia_path = "/home/energy/tipapa/Yang_work/Pt/"
coordination_sphere_1 = ["3", "4"]  # everything else has weird pores in them
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
        co_str = f"N{sphere_one_ount}+{sphere_two_count}N+{sphere_three_count}N"
        n_count = sphere_one_ount + sphere_two_count + sphere_three_count
        n_counter.append(n_count)
        print(co_str, n_count)

    print(Counter(n_counter))


# look_at_data()

def check_coordination_sphere(structure):
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
    return sphere_one_ount, sphere_two_count, sphere_three_count

# def get_structures():
#     for file in files:
#         if "png" in file:
#             dir = file.split("/")[-1].replace(".png", "")
#             try:
#                 os.mkdir(dir)
#             except:
#                 pass
#             structure_path = (
#                 Path(tia_path) / file.split("/")[-1].replace(".png", "") / "POSCAR.opt"
#             )
#             os.system(f"cp {structure_path} {dir}")
#             structure = read(structure_path)
#             write(f"./{dir}/{file.split('/')[-1]}", structure)

# get_structures()

def replace_first_digit(input_string, replacement_letter, count):
    # Use a function to replace the first match
    def replacer(match):
        return replacement_letter

    # Replace the first digit using re.sub with a count of 1
    result = re.sub(r'\d', replacer, input_string, count=count)
    return result

def make_dirs(file, structures, coor_sphere):
    nsp1, nsp2, nsp3 = check_coordination_sphere(structures)
    sp1, sp2, sp3 = check_coordination_sphere(coor_sphere)
    # print(nsp1, nsp2, nsp3)
    # print(sp1, sp2, sp3)
    if nsp1 != sp1:
        b = nsp1 - sp1
        filename = file.split("/")[-2].replace(".png", "")
        new_string = f'{nsp1-b}B{b}'
        output_string  = replace_first_digit(filename, new_string, count=1)
        if '+' in filename:
            fileend = filename.split('+')[1]
            try:
                os.mkdir(f'{parent_dir}{output_string}')
                write(f'{parent_dir}{output_string}/POSCAR.opt', coor_sphere)
                write(f'{parent_dir}{output_string}/{output_string}.png', coor_sphere)
            except:
                pass
        else:
            try:
                os.mkdir(f'{parent_dir}{output_string}')
                write(f'{parent_dir}{output_string}/POSCAR.opt', coor_sphere)
                write(f'{parent_dir}{output_string}/{output_string}.png', coor_sphere)
            except:
                pass

    if '+' in file:          
        if nsp2 != sp2 or nsp3 != sp3:
            diff1 = nsp2 - sp2
            diff2 = nsp3 - sp3            
            n_content = nsp2 + nsp3
            if diff1 + diff2 == 1:
                b = 1
            else:
                b = 2
            new_string = f'{n_content-b}N{b}B'
        else:
            new_string = file.split('/')[-2].split('+')[1]
        if nsp1 == sp1:
            filestart = filename.split('+')[0]
        else:
            bb = nsp1 - sp1
            first_string = f'{nsp1-bb}B{bb}'
            filestart  = replace_first_digit(filename, first_string, count=1)
            
        filename = file.split("/")[-2].replace(".png", "")
        if '-' in filename:
            fileend = filename.split('-')[1]
            ss = filestart.split('+')[0]
            string = f'{ss}+{new_string}-{fileend}'
            try:
                os.mkdir(f'{parent_dir}{string}')
            except:
                pass
            write(f'{parent_dir}{string}/POSCAR.opt', coor_sphere)
            poscar = read(f'{parent_dir}{string}/POSCAR.opt')
            write(f'{parent_dir}{string}/{string}.png', poscar)
        else:
            ss = filestart.split('+')[0]
            string = f'{ss}+{new_string}'
            try:
                os.mkdir(f'{parent_dir}{string}')
            except:
                pass
            write(f'{parent_dir}{string}/POSCAR.opt', coor_sphere)
            poscar = read(f'{parent_dir}{string}/POSCAR.opt')
            write(f'{parent_dir}{string}/{string}.png', poscar)
    print(f'{parent_dir}{string}')

dictss = {}
def create_doped_structs():
    for file in files[:]:
        if '+' in file:
            print(file)
            n_idx = []
            structures = read(file)
            for idx, struc in enumerate(structures.get_chemical_symbols()):
                if struc == "N":
                    n_idx.append(idx)        
            view(structures)
            inp1 = int(input("enter indices only 1 indice  ")) # indice for structure with 1 boron sub
            inp2 = input("enter indices only 2 indice  ")  # indice for structure with 2 boron subs
            # os.system(f'ase gui {file}')
            # os.system('echo ronan')
            b1 = int(inp2.split(" ")[0])
            # b1 = 35
            # b2 = 34
            # inp1 = 35
            b2 = int(inp2.split(" ")[1])
            indice1_strucs = structures.copy()
            indice2_strucs = structures.copy()

            if inp1:
                indice1_strucs[inp1].symbol = "B"
                print(indice1_strucs[inp1].symbol)
            if b1:
                indice2_strucs[b1].symbol = "B"
                indice2_strucs[b2].symbol = "B"

            make_dirs(file, structures, indice1_strucs)
            make_dirs(file, structures, indice2_strucs)

create_doped_structs()
