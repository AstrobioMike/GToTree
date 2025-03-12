from setuptools import setup, find_packages

setup(
    name="GToTree",
    version="2.0.0",
    description="A user-friendly workflow for phylogenomics",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="AstrobioMike",
    url="https://github.com/AstrobioMike/GToTree",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "gtotree2=gtotree.cli.main:main",
        ],
    },
    install_requires=[
        "biopython",
        "pyarrow",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GPL3",
    ],
    python_requires='==3.12.7',
)
