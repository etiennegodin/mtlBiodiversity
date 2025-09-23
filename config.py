import os
from dataclasses import dataclass
from pathlib import Path
from dotenv import load_dotenv


load_dotenv()  # reads .env if present
@dataclass
class Config:
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


configs = load_config()