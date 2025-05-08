from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_desc = fh.read()
setup(
    name="GToTree",
    version="2.0.0",
    description="A user-friendly workflow for phylogenomics",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="AstrobioMike",
    url="https://github.com/AstrobioMike/GToTree",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "gtotree2=gtotree.main:main",
            "gtt-midpoint-root-tree=gtotree.utils.helper_scripts.gtt_midpoint_root_tree:main",
        ],
    },
    install_requires=[
        "biopython",
        "pyarrow",
        "pyhmmer",
        "tqdm",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GPL3",
    ],
    python_requires='==3.12.7',
)
