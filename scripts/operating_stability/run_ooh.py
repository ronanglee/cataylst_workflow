# """Code taken and refactored from Tiaporn"""

# import glob
# import os


# def get_wd(folders):
#     cwd = os.getcwd()
#     for folder in folders:
#         cwd = os.path.join(cwd, folder)
#     return cwd


# def electronic_done(cwd):
#     finish = 0
#     for item in os.listdir(cwd):
#         if item == "ACF.dat":
#             finish = 1
#             break
#     return finish


# home_dir = os.getcwd()


# def main(**data: dict) -> None:

#     xc_list = ["BEEF-vdW"]

#     prototype_list = get_all_sacs_stability()
#     dopant_list = ["non"]
#     dopantsite_list = ["non"]
#     solvation_list = ["implicit", "vac"]
#     ads2_list = ["non", "OOH"]
#     for xc in xc_list:
#         for prototype in prototype_list:
#             for metal in metal_list:
#                 prototype = prototype.replace("Pt", "M")
#                 for dopant in dopant_list:
#                     for dopantsite in dopantsite_list:
#                         ads1_list = ["non"]
#                         for ads1 in ads1_list:
#                             if ads1 == "non":
#                                 ads1site = "non"
#                                 ads1orient = "non"
#                             for ads2 in ads2_list:
#                                 if ads2 == "non":
#                                     ads2site = "non"
#                                     ads2orient = "non"
#                                 else:
#                                     ads2site = "ontop_metal"
#                                     ads2orient = "1"
#                                 for solvation in solvation_list:
#                                     if solvation == "implicit" and ads2 != "non":
#                                         run_folders = [
#                                             "vibration",
#                                             "vasp_rx",
#                                         ]  # or vibration
#                                     else:
#                                         run_folders = ["vasp_rx"]
#                                     for run_folder in run_folders:
#                                         user_args = [
#                                             xc,
#                                             prototype,
#                                             metal,
#                                             "1",
#                                             dopant,
#                                             dopantsite,
#                                             ads1,
#                                             ads1site,
#                                             ads1orient,
#                                             ads2,
#                                             ads2site,
#                                             ads2orient,
#                                             solvation,
#                                             run_folder,
#                                         ]
#                                         default = user_args
#                                         folder = get_wd(default)
#                                         if (
#                                             run_folder == "vibration"
#                                             and ads1 == "non"
#                                             and ads2 == "non"
#                                         ):
#                                             continue
#                                         os.makedirs(folder, exist_ok=True)
#                                         os.chdir(folder)
#                                         cwd = os.getcwd()
#                                         # find = electronic_done(cwd)
#                                         # if find == 1:
#                                         #    print('Done')
#                                         #    from ase.io import read, write
#                                         os.system(
#                                 "python /home/energy/rogle/asm_orr_rxn/generated_data/tia_MN4/build_initial.py"
#                                         )
#                                         #    os.system('python /home/energy/tipapa/tools/extract_features_graph.py')
#                                         # else:
#                                         #    print('Error --- Not Done')
#                                         os.chdir(home_dir)
