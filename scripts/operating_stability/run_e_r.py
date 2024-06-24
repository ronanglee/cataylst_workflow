# """Code taken and refactored from Tiaporn"""

# import glob
# import json
# import os
# import re

# import matplotlib.pyplot as plt
# import numpy as np
# from ase.db import connect
# from ase.io import read
# from get_ooh import main_ooh
# from matplotlib.patches import Patch
# from zero_h_synthesis_vals import get_gradient_vals


# def call_db(sac):
#     db = connect(f"/home/energy/rogle/datasets/orr_complete_data/oh_databases/{sac}.db")
#     return db


# with open(
#     "/home/energy/rogle/asm_orr_rxn/check_if_structures_are_avoided/avoid.json", "r"
# ) as f:
#     avoid_motifs = json.load(f)

# with open(r"ooh_etc_thermal_corrections.json", "r") as f:
#     thermal_corrections = json.load(f)

# with open(r"new_corrected_carbon_sheet_thermal_corrections.json", "r") as f:
#     corrected_carbon_sheet_thermal_corrections = json.load(f)

# # constant
# kB = 8.617e-5  # Boltzmann constant in eV/K
# h = 4.136e-15  # Planck constant in eV*s
# T = 298.15  # Temperature in K
# k = kB * T * np.log(10)
# conc = 10**-6  # concentration of metal ion
# RT = 0.0257149  # at T = 298.15K

# # gas molecules at T = 298.15 K, p = 1 bar
# # eV/molecule (EDFT + ZPE + H - TS + BEEF correction)
# G_H2 = -7.128  # eV/molecule (Free energy + BEEF correction)
# G_H2O = -12.869  # eV/molecule (Free energy + BEEF correction)


# def acid_stability(
#     pH,
#     u,
#     m,
#     G_MN4,
#     G_CH0,
#     G_CH1,
#     G_CH2,
#     G_CH3,
#     G_CH4,
#     G_MN4_H2O,
#     G_MN4_OH,
#     G_MN4_O,
#     G_MN4_H,
#     G_MN4_OOH,
# ):
#     if m == "Pd":
#         G_M = -1.96  # ev/atom
#         G_M2 = G_M + 2 * (0.951 + (1 / 2) * 0.0592 * np.log10(conc))  # 1
#         # 1. dissolution to M+2
#         dg1 = G_M2 + G_CH0 - G_MN4 - 2 * u - 0.0 * G_H2 + 0 * k * pH
#         dg2 = G_M2 + G_CH1 - G_MN4 - 1 * u - 0.5 * G_H2 + 1 * k * pH
#         dg3 = G_M2 + G_CH2 - G_MN4 - 0 * u - 1.0 * G_H2 + 2 * k * pH
#         dg4 = G_M2 + G_CH3 - G_MN4 + 1 * u - 1.5 * G_H2 + 3 * k * pH
#         dg5 = G_M2 + G_CH4 - G_MN4 + 2 * u - 2.0 * G_H2 + 4 * k * pH
#         dg6 = dg7 = dg8 = dg9 = dg10 = 9999
#         dg11 = dg12 = dg13 = dg14 = dg15 = 9999
#         dg16 = dg17 = dg18 = dg19 = dg20 = 9999
#         dg21 = dg22 = dg23 = dg24 = dg25 = 9999
#         dg26 = dg27 = dg28 = dg29 = dg30 = 9999
#         dg31 = dg32 = dg33 = dg34 = dg35 = 9999
#         dg36 = dg37 = dg38 = dg39 = dg40 = 9999
#     elif m == "Pt":
#         G_M = -3.13  # eV/atom
#         G_M2 = G_M + 2 * (1.18 + (1 / 2) * 0.0592 * np.log10(conc))  # 1
#         # 1. dissolution to M+2
#         dg1 = G_M2 + G_CH0 - G_MN4 - 2 * u - 0.0 * G_H2 + 0 * k * pH
#         dg2 = G_M2 + G_CH1 - G_MN4 - 1 * u - 0.5 * G_H2 + 1 * k * pH
#         dg3 = G_M2 + G_CH2 - G_MN4 - 0 * u - 1.0 * G_H2 + 2 * k * pH
#         dg4 = G_M2 + G_CH3 - G_MN4 + 1 * u - 1.5 * G_H2 + 3 * k * pH
#         dg5 = G_M2 + G_CH4 - G_MN4 + 2 * u - 2.0 * G_H2 + 4 * k * pH
#         dg6 = dg7 = dg8 = dg9 = dg10 = 9999
#         dg11 = dg12 = dg13 = dg14 = dg15 = 9999
#         dg16 = dg17 = dg18 = dg19 = dg20 = 9999
#         dg21 = dg22 = dg23 = dg24 = dg25 = 9999
#         dg26 = dg27 = dg28 = dg29 = dg30 = 9999
#         dg31 = dg32 = dg33 = dg34 = dg35 = 9999
#         dg36 = dg37 = dg38 = dg39 = dg40 = 9999

#     dg0 = 0
#     dg41 = G_MN4_H2O - G_MN4 - G_H2O
#     dg42 = G_MN4_OH + 0.5 * G_H2 - 1 * k * pH - u - G_MN4 - G_H2O
#     dg43 = G_MN4_O + G_H2 - 2 * k * pH - 2 * u - G_MN4 - G_H2O
#     dg44 = G_MN4_H - G_MN4 - 0.5 * G_H2 + 1 * k * pH + u
#     # ooh + (1.5 * G_H2) - np - (2*G_H2O) + data['data'][string] + 0.2
#     dg45 = G_MN4_OOH + (1.5 * G_H2) - 1 * k * pH - 1 * u - G_MN4 - (2 * G_H2O)

#     dg = [
#         dg0,
#         dg1,
#         dg2,
#         dg3,
#         dg4,
#         dg5,
#         dg6,
#         dg7,
#         dg8,
#         dg9,
#         dg10,
#         dg11,
#         dg12,
#         dg13,
#         dg14,
#         dg15,
#         dg16,
#         dg17,
#         dg18,
#         dg19,
#         dg20,
#         dg21,
#         dg22,
#         dg23,
#         dg24,
#         dg25,
#         dg26,
#         dg27,
#         dg28,
#         dg29,
#         dg30,
#         dg31,
#         dg32,
#         dg33,
#         dg34,
#         dg35,
#         dg36,
#         dg37,
#         dg38,
#         dg39,
#         dg40,
#         dg41,
#         dg42,
#         dg43,
#         dg44,
#         dg45,
#     ]
#     return dg


# def acid_stability_H2O(
#     pH,
#     u,
#     m,
#     G_MN4,
#     G_MN4_O,
#     G_MN4_OH,
#     G_CH0,
#     G_CH1,
#     G_CH2,
#     G_CH3,
#     G_CH4,
#     G_MN4_H2O,
#     G_MN4_H,
#     G_MN4_OOH,
# ):
#     dg = acid_stability(
#         pH,
#         u,
#         m,
#         G_MN4,
#         G_CH0,
#         G_CH1,
#         G_CH2,
#         G_CH3,
#         G_CH4,
#         G_MN4_H2O,
#         G_MN4_OH,
#         G_MN4_O,
#         G_MN4_H,
#         G_MN4_OOH,
#     )
#     react_dict = {
#         dg[1]: 1,
#         dg[2]: 1,
#         dg[3]: 1,
#         dg[4]: 1,
#         dg[5]: 1,  # 1
#         dg[6]: 2,
#         dg[7]: 2,
#         dg[8]: 2,
#         dg[9]: 2,
#         dg[10]: 2,  # 2
#         dg[11]: 3,
#         dg[12]: 3,
#         dg[13]: 3,
#         dg[14]: 3,
#         dg[15]: 3,  # 3
#         dg[16]: 4,
#         dg[17]: 4,
#         dg[18]: 4,
#         dg[19]: 4,
#         dg[20]: 4,  # 4
#         dg[21]: 5,
#         dg[22]: 5,
#         dg[23]: 5,
#         dg[24]: 5,
#         dg[25]: 5,  # 5
#         dg[26]: 6,
#         dg[27]: 6,
#         dg[28]: 6,
#         dg[29]: 6,
#         dg[30]: 6,  # 6
#         dg[31]: 7,
#         dg[32]: 7,
#         dg[33]: 7,
#         dg[34]: 7,
#         dg[35]: 7,  # 7
#         dg[36]: 8,
#         dg[37]: 8,
#         dg[38]: 8,
#         dg[39]: 8,
#         dg[40]: 8,  # 8
#         dg[0]: 9,  # bare
#         dg[41]: 10,  # H2O
#         dg[42]: 11,  # OH
#         dg[43]: 12,  # O
#         dg[44]: 13,  # H
#         dg[45]: 14,
#     }  # OOH

#     react = react_dict[min(react_dict)]
#     return react


# def relative_stability_H2O(
#     pH,
#     u,
#     m,
#     G_MN4,
#     G_MN4_O,
#     G_MN4_OH,
#     G_CH0,
#     G_CH1,
#     G_CH2,
#     G_CH3,
#     G_CH4,
#     G_MN4_H2O,
#     G_MN4_H,
#     G_MN4_OOH,
# ):
#     dg = acid_stability(
#         pH,
#         u,
#         m,
#         G_MN4,
#         G_CH0,
#         G_CH1,
#         G_CH2,
#         G_CH3,
#         G_CH4,
#         G_MN4_H2O,
#         G_MN4_OH,
#         G_MN4_O,
#         G_MN4_H,
#         G_MN4_OOH,
#     )
#     # relative energy
#     ref = min(dg[1:41])
#     B_rel = dg[0] - ref
#     H2O_rel = dg[41] - ref
#     OH_rel = dg[42] - ref
#     O_rel = dg[43] - ref
#     H_rel = dg[44] - ref
#     OOH_rel = dg[45] - ref
#     react_dict = {
#         B_rel: 0,  # bare
#         #   H2O_rel: 1,  # *H2O ## not doing water
#         OH_rel: 2,  # *OH
#         O_rel: 3,  # *O
#         H_rel: 4,
#         OOH_rel: 5,
#     }  # *H
#     react_border = react_dict[min(react_dict)]
#     react = min(B_rel, H2O_rel, OH_rel, O_rel, H_rel)
#     return (react, react_border, B_rel, H2O_rel, OH_rel, O_rel, H_rel, OOH_rel)


# def get_data_corrected_carbon(dictin, sac, label):
#     sac = sac.replace("Pt", "M")
#     if label == "CH0":
#         label = None
#         joining = f"{sac}"
#     else:
#         label = f"_{label}"
#         joining = f"{sac}{label}"
#     return dictin["data"][joining]


# def get_data_thermal_corrections(dictin, sac, metal, label1, label2):
#     sac = sac.replace("Pt", "M")
#     joining = f"{sac}_{metal}_{label1}_{label2}"
#     return dictin["data"][joining]


# def get_coordination_spheres():
#     sacs = glob.glob(
#         f"/home/energy/rogle/asm_orr_rxn/local_structure/nanocluster_formation/data/e_mnxc/Pt/PtN*"
#     )
#     used_sacs = glob.glob(
#         f"/home/energy/rogle/asm_orr_rxn/local_structure/nanocluster_formation/data/e_nch_metal_removed/*"
#     )  # master directory
#     base_sacs = [i.split("/")[-1] for i in used_sacs]

#     coordination_strings = []
#     for sac in sacs:
#         if sac.split("/")[-1] not in base_sacs:
#             continue
#         structure = read(f"{sac}/POSCAR.opt")
#         m_idx = [i for i, x in enumerate(structure.get_chemical_symbols()) if x == "Pt"]
#         n_idx = [i for i, x in enumerate(structure.get_chemical_symbols()) if x == "N"]
#         distances = []
#         for i in m_idx:
#             for j in n_idx:
#                 distances.append(structure.get_distance(i, j))

#         sphere_one_ount = 0
#         sphere_two_count = 0
#         sphere_three_count = 0
#         for i in sorted(distances):
#             if i < 2.6:
#                 sphere_one_ount += 1
#             elif i < 3.6 and i > 2.6:
#                 sphere_two_count += 1
#             elif i > 3.6:
#                 sphere_three_count += 1
#         co_str = f"PtN{sphere_one_ount}+{sphere_two_count}N+{sphere_three_count}N"
#         coordination_strings.append(co_str)

#     return coordination_strings


# def get_all_sacs_stability():
#     path = "/home/energy/rogle/asm_orr_rxn/local_structure/nanocluster_formation/data/e_nch_metal_removed"
#     sacs_path = glob.glob(f"{path}/*")
#     sacs = []
#     for sac in sacs_path:
#         basename = os.path.basename(sac)
#         basename = basename.replace("Pt", "M")
#         sacs.append(basename)
#     return sacs


# def get_gradient_array(carbon, m):
#     gradient_array = []
#     for C in carbon:
#         gradient_array.append(get_gradient_vals(C, m))

#     x_max = max(gradient_array, key=lambda x: x[0])
#     y_max = max(gradient_array, key=lambda x: x[1])
#     x_min = min(gradient_array, key=lambda x: x[0])
#     y_min = min(gradient_array, key=lambda x: x[1])

#     euc_dis_array = []
#     dis_max_min = np.sqrt((x_max[0] - x_min[0]) ** 2 + (y_max[1] - y_min[1]) ** 2)
#     for values in gradient_array:
#         # find euchlidean distance
#         euc_dis_array.append(
#             np.sqrt((values[0] - x_max[0]) ** 2 + (values[1] - y_max[1]) ** 2)
#             / dis_max_min
#         )
#     return euc_dis_array


# carbon = get_all_sacs_stability()
# coordination_spheres = get_coordination_spheres()

# electrolyte = "Bare"

# u_orr = 0.7
# pH_orr = 0

# colors = ["b", "g", "r", "c", "m", "black", "y", "green"]
# shapes = [
#     ".",
#     "o",
#     "v",
#     "^",
#     "<",
#     ">",
#     "1",
#     "2",
#     "3",
#     "4",
#     "s",
#     "p",
#     "*",
#     "h",
#     "H",
#     "+",
#     "x",
#     "D",
#     "d",
#     "|",
#     "_",
#     "P",
#     "X",
#     "8",
#     "d",
#     "v",
#     "<",
#     ">",
#     "^",
# ]
# nitrogen_dict = {3: "blue", 4: "green", 5: "red", 6: "orange", 7: "purple", 8: "brown"}

# metals = ["Pd", "Pt"]

# color_dict = {}
# shapes_dict = {}
# for col, m in zip(colors, metals):
#     color_dict[m] = col

# for C, sphere, shape in zip(carbon, coordination_spheres, shapes):
#     Gr_orr = []
#     Gr_orr_name = []
#     ooh_vals = []
#     db = call_db(C)
#     numbers = re.findall(r"\d+", sphere)
#     n = sum([int(i) for i in numbers])
#     if "A" in C:
#         extra_string = f"{m}N{n} (A - {sphere.replace('Pt', '')})"
#     elif "Z" in C:
#         extra_string = f"{m}N{n} (Z - {sphere.replace('Pt', '')})"
#     else:
#         extra_string = f"{m}N{n} ({sphere.replace('Pt', '')})"
#     for key in shapes_dict.keys():
#         if extra_string == key:
#             indices = len(
#                 [
#                     index
#                     for index, string in enumerate(shapes_dict.keys())
#                     if extra_string in string
#                 ]
#             )
#             extra_string += f" - {indices+1}"
#     shapes_dict[extra_string] = shape


# def main(metal):
#     fig, ax = plt.subplots()
#     annote = []
#     euc_dis_array = get_gradient_array(carbon, metal)
#     for idx, C in enumerate(carbon):
#         ooh_vals = []
#         db = call_db(C)
#         continue_loop = False
#         for key, vals in avoid_motifs["data"].items():
#             m_c = (
#                 "MN" + C.split("N", 1)[1]
#             )  # do only max split of 1 for first occurance of N
#             if key == m_c and vals == metal:
#                 print(m_c, metal, "avoided")
#                 continue_loop = True
#         if continue_loop:
#             continue
#         G_CH0 = get_data_corrected_carbon(
#             corrected_carbon_sheet_thermal_corrections, C, "CH0"
#         )
#         G_CH1 = get_data_corrected_carbon(
#             corrected_carbon_sheet_thermal_corrections, C, "1H"
#         )
#         try:  # MN3CZ_2H
#             G_CH2 = get_data_corrected_carbon(
#                 corrected_carbon_sheet_thermal_corrections, C, "2H"
#             )
#         except:
#             G_CH2 = 9999
#         try:  # MN6C-2+1N-2_3H not present
#             G_CH3 = get_data_corrected_carbon(
#                 corrected_carbon_sheet_thermal_corrections, C, "3H"
#             )
#         except:
#             G_CH3 = 9999
#         try:
#             G_CH4 = get_data_corrected_carbon(
#                 corrected_carbon_sheet_thermal_corrections, C, "4H"
#             )
#         except:  # not possible to get CH4
#             G_CH4 = 9999
#         thermal_OOH = get_data_thermal_corrections(
#             thermal_corrections, C, metal, "non", "OOH"
#         )
#         thermal_H = 9999
#         gr_orr = []
#         gr_orr_name = []
#         try:
#             E_MN4 = db.get(metal=metal, ads1="non", ads2="non").energy
#         except:
#             print(f"carbon:{C} metal:{metal} not present")
#             continue

#         # *OH
#         E_MN4_OH = 9999

#         # *OOH

#         E_MN4_OOH = db.get(metal=metal, ads1="non", ads2="OOH").energy

#         # *O
#         E_MN4_O = 9999

#         G_MN4 = E_MN4
#         # G_MN4_H2O = E_MN4_H2O + thermal_H2O
#         G_MN4_H2O = 9999
#         G_MN4_OH = 9999
#         G_MN4_O = 9999
#         G_MN4_H = 9999
#         G_MN4_OOH = E_MN4_OOH + thermal_OOH + 0.2

#         temp1, temp1_border, B_rel, H2O_rel, OH_rel, O_rel, H_rel, OOH_rel = (
#             relative_stability_H2O(
#                 pH_orr,
#                 u_orr,
#                 m,
#                 G_MN4,
#                 G_MN4_O,
#                 G_MN4_OH,
#                 G_CH0,
#                 G_CH1,
#                 G_CH2,
#                 G_CH3,
#                 G_CH4,
#                 G_MN4_H2O,
#                 G_MN4_H,
#                 G_MN4_OOH,
#             )
#         )
#         if electrolyte == "Bare":
#             electro = B_rel
#         elif electrolyte == "OOH":
#             electro = OOH_rel
#         gr_orr.append(electro)
#         gr_orr_name.append("*")
#         vals = main_ooh(C, metal)
#         ooh_vals.append(vals)
#         # if vals > 3.7 and vals < 4.7:
#         #     annote.append(ax.text(vals, gr_orr[0], f'{m}'))
#         n_str = list(shapes_dict.keys())[idx].split("(", 1)[1].split(")")[0]
#         numbers = re.findall(r"\d+", n_str)
#         n = sum([int(i) for i in numbers])
#         p1 = plt.scatter(
#             vals,
#             gr_orr,
#             c=nitrogen_dict[n],
#             marker=list(shapes_dict.values())[idx],
#             s=75,
#         )

#     def f(c, m):
#         return plt.plot([], [], marker=m, color=c, ls="none")[0]

#     handles = []
#     n_used = []
#     for key, val in shapes_dict.items():
#         n_str = key.split("(", 1)[1].split(")")[0]
#         # print(n_str)
#         numbers = re.findall(r"\d+", n_str)
#         n = sum([int(i) for i in numbers])
#         if n not in n_used:
#             n_used.append(n)
#         handles.append(f(nitrogen_dict[n], val))

#     legend_patches = []

#     for nitro in sorted(n_used):
#         legend_patches.append(Patch(color=nitrogen_dict[nitro], label=f"N={nitro}"))

#     ax.legend(handles, list(shapes_dict.keys()), bbox_to_anchor=(1.05, 1), ncol=2)
#     first_legend = ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 0.3))
#     fig.gca().add_artist(first_legend)
#     # add second legend

#     legend = plt.legend(
#         handles, list(shapes_dict.keys()), bbox_to_anchor=(1.05, 1), ncol=2
#     )
#     # legend_patches = [Patch(color=color, label=name) for name, color in color_dict.items()]
#     legend.set_title("# atoms (configuration)")
#     ax.axvline(x=4.2, color="r", linestyle="--", linewidth=2)
#     # first_legend = ax.legend(handles=legend_patches, bbox_to_anchor=(1.2, 0.6))
#     # ax.legend(bbox_to_anchor=(1.2, 0.9))
#     # ax.legend(loc='lower left', ncol=2, bbox_to_anchor=(1.15, 0.3))
#     # ax.add_artist(first_legend)
#     ax.set_ylabel("$\Delta$G$_R$ \ eV")
#     plt.title(f"{electrolyte} surface site @{u_orr} V$_{{SHE}}$, 0 pH")

#     ax.set_xlabel("$\Delta$G(OOH*) \ eV (theoretical)")
#     plt.savefig(
#         f"nitrogens_{metal}_zero_h_gradient_{electrolyte}_stability_plot_{u_orr}.png",
#         bbox_inches="tight",
#     )


# for met in metals:
#     main(met)
