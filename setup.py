from setuptools import find_packages, setup

setup(
    name="phylotres",
    version="0.0.1",
    keywords=("pip", "phylotres"),
    description="PhyloTres",
    long_description="PhyloTres",
    license="GNU GENERAL V3.0",
    url="https://github.com/cribbslab/phylotres; https://www.ndorms.ox.ac.uk/research/research-groups/cribbs-group-computational-and-systems-biology",
    author="Jianfeng Sun",
    author_email="jianfeng.sun@ndorms.ox.ac.uk; Adam.Cribbs@ndorms.ox.ac.uk",
    packages=find_packages(),
    include_package_data=True,
    # package_data={},
    platforms="any",
    python_requires=">=3.10",
    install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "biopython",
        "pyfiglet==0.8.post1",
    ],
)