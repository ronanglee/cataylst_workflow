# import os
# from pathlib import Path

# from ase.io import read  # type: ignore


# def main(base_dir: os.PathLike) -> None:
#     """Calculate the operating stability of the catalyst.
#     g_a = e_xc + e_m_on_c - e_c - e_mxc
#     g_d = e_m + e_xc - e_mxc
#     where x is the dopant, c is carbon and m is the metal.

#     Args:
#         base_dir (Path): Path to the base directory.

#     Returns:
#         None
#     """
#     data_base_folder = Path(base_dir) / "runs" / "databases"
#     e_m_on_c_db = data_base_folder / "e_m_on_c.db"
#     e_c_db = data_base_folder / "e_c.db"
#     e_mxc_db = data_base_folder / "e_mxc.db"
#     e_m_db = data_base_folder / "e_m.db"
#     structures = []
#     for row in e_mxc_db.select():
#         structures.append(row.name)
