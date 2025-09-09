import pandas as pd 

df = pd.read_parquet("data/processed/annual_observations.parquet")

print(df)