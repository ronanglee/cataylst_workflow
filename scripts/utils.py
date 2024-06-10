import os
from pathlib import Path

from ase.db import connect  # type: ignore
from ase.io import read  # type: ignore


def read_and_write_datbase(struc_dir: os.PathLike, base_dir: os.PathLike) -> None:
    """Read and write the database.

    Args:
        struc_dir (Path): Path to the directory where to find the initial geometry.
        base_dir (Path): Path to the base directory.

    Returns:
        None

    """
    data_base_folder = Path(base_dir) / "runs" / "databases"
    name = data_base_folder.stem
    db = connect(data_base_folder / f"{name}.db")
    structure = read(os.path.join(struc_dir, "OUTCAR"))
    db.write(structure, key_value_pairs={"name": name})
