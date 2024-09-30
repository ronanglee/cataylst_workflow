import glob
import os
import sqlite3
from pathlib import Path  # type: ignore

from ase.db import connect  # type: ignore


def insert_data(save_file: os.PathLike, data: tuple) -> None:
    """Insert data into the database.

    Args:
        save_file (os.PathLike): Path to the database.
        data (tuple): Data to insert into the database.

    """
    conn = sqlite3.connect(f"{save_file}")
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

    tables = cursor.fetchall()

    for table in tables:
        name = table[0]
        break
    cursor.execute(f"PRAGMA table_info({name})")
    columns = cursor.fetchall()
    # columns[1] gets the second column, [1] gets the name of the column
    second_column_name = columns[1][1]
    check_query = (
        f"SELECT 1 FROM {name} WHERE {second_column_name} = ? AND value = ? LIMIT 1"
    )
    if isinstance(data[1], bytes):
        check_data = (data[0], data[1].decode("utf-8"))
    cursor.execute(check_query, (data[0], data[1]))

    result = cursor.fetchone()
    if result:
        pass
    else:
        print(f"inserting into {save_file}", check_data, flush=True)
        config_name, value = check_data
        cursor.execute(
            """
            INSERT INTO Configurations (config_name, value)
            VALUES (?, ?)
        """,
            (config_name, value),
        )

    conn.commit()
    conn.close()


def read_data(file: os.PathLike) -> list:
    """Read data from the database.

    Args:
        file (os.PathLike): Path to the database.

    Returns:
        data (list): List of data from the database.
    """
    conn = sqlite3.connect(file)
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    for table in tables:
        name = table[0]
        break

    cursor.execute(f"SELECT * FROM {name}")
    data = cursor.fetchall()
    conn.close()
    return data


def merge(local_db: os.PathLike, master_db: os.PathLike) -> None:
    """Merge the local database with the master database.

    Args:
        local_db (os.PathLike): Path to the local database.
        master_db (os.PathLike): Path to the master database.

    """
    dbs = glob.glob(os.path.join(master_db, "*.db"))
    print("reading", dbs, flush=True)
    for db in dbs:
        if "e_m_" in db:
            continue
        print(f"checking {db} in {local_db}", flush=True)
        if "vib" in db and "corrections" in db:
            indiv_db = os.path.join(local_db, db.replace("_master", ""))
            if not os.path.exists(indiv_db):
                print(f"{indiv_db} does not exist", flush=True)
                continue
            data = read_data(Path(indiv_db))
            for d in data:
                if isinstance(d[2], float):
                    insert_data(Path(db), (d[1], d[2]))
                else:
                    key1 = d[1]
                    key2 = d[2]
                    insert_data(Path(db), (key1, key2))
        else:
            new_loc = os.path.join(local_db, db.replace("_master", ""))
            if os.path.exists(new_loc):
                master_db_connected = connect(db)
                local_db_connected = connect(new_loc)
                for row in local_db_connected.select():
                    try:
                        master_db_connected.get(name=row.name)
                    except NameError:
                        print(f"inserting {row.name} into {db}", flush=True)
                        master_db_connected.write(row)
