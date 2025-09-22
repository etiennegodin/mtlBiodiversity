import os
from dataclasses import dataclass
from dotenv import load_dotenv
from mtlBio.core import raw_data_read_only
from pathlib import Path


@dataclass
class Config:
    raw_data_read_only('data/raw', debug = False)
    data_dir: str
    sql_dir : str
    limit : int
    debug : bool

def load_config() -> Config:
    return Config(
        data_dir = Path(os.getenv("DATA_PATH")),
        sql_dir = Path(os.getenv("SQL_PATH")),
        limit = os.getenv('LIMIT'),
        debug = os.getenv('DEBUG')
    )

# Single shared instance
configs = load_config()