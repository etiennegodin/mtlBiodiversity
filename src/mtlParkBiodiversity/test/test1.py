import pandas as pd 
import duckdb
from matplotlib import pyplot as plt
import seaborn as sns
df = pd.read_parquet("data/processed/most_observed_species.parquet")

df = df.set_index('species')


