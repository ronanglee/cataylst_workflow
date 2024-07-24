import glob
from collections import Counter
from ase.io import read


print(len(glob.glob('test_structs/Pt*')))
print(len(glob.glob('test_structs/dopant/Pt*')))
# n_content = []
# b_content = []
# for f in glob.glob('test_structs/dopant/*/POSCAR.opt'):
#     print(f)
#     atoms = read(f)
#     counter = dict(Counter(atoms.get_chemical_symbols()))
#     print(counter)
#     print(len(atoms))
#     n_content.append(round(counter['N'], 0))
#     b_content.append(round(counter['B'], 0))

# print(set(n_content))
# print(set(b_content))
# print(Counter(n_content))
# print(Counter(b_content))