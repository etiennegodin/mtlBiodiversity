from setuptools import setup, find_packages

setup(
    name='mtlParkBiodiversity',
    version='0.0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'geopandas',
        'pandas',
        'numpy',
        'duckdb',
        'duckdb-extension-spatial',
        'pyarrow'
    ],
    )