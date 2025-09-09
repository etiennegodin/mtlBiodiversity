from setuptools import setup, find_packages

setup(
    name='mtlParkBiodiversity',
    version='0.0.2',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points = {
        "console_scripts": "mbio = mtlParkBiodiversity.cli:main"
    },
    install_requires=[
        'geopandas',
        'pandas',
        'numpy',
        'duckdb',
        'duckdb-extension-spatial',
        'fiona'
    ],
    )