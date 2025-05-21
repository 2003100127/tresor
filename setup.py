from setuptools import find_packages, setup

setup(
    name="tresor",
    version="0.1.1",
    keywords=["conda", "tresor"],
    description="Tresor",
    long_description="Simulation tool Tresor",
    license="GNU-3.0",

    url="https://github.com/2003100127/tresor",
    author="Jianfeng Sun",
    author_email="jianfeng.sun@ndorms.ox.ac.uk",

    packages=find_packages(),
    include_package_data=True,
    # package_data={},
    platforms="any",
    python_requires=">3.9",
    install_requires=[
        "click",
        "pandas",
        "numpy",
        "scipy",
        "biopython",
        "tables",
        "pyyaml",
        # "rpy2", # 3.4.5
        "pyfiglet==0.8.post1",
    ],
    entry_points={
        'console_scripts': [
            'tresor=tresor.main:main',
        ],
    },
)