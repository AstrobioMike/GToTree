from setuptools import setup

setup(
    name="GToTree",
    version="1.4.2",
    url="https://github.com/AstrobioMike/GToTree",
    project_urls={
        "Documentation": "https://github.com/AstrobioMike/GToTree/wiki",
        "Code": "https://github.com/AstrobioMike/GToTree",
    },
    license="GNU General Public License v3.0",
    author="Michael D. Lee",
    author_email="Mike.Lee@nasa.gov",
    description="A user-friendly workflow for phylogenomics.",
    long_description="GToTree is a user-friendly workflow for phylogenomics intended to give more researchers the capability to create phylogenomic trees. The open-access Bioinformatics Journal publication is available here (https://doi.org/10.1093/bioinformatics/btz188), and documentation and examples can be found at the wiki here (https://github.com/AstrobioMike/GToTree/wiki).",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
    scripts=['bin/*'],
)
