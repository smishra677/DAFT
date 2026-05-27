from pathlib import Path
from setuptools import find_namespace_packages, setup

ROOT = Path(__file__).parent
README = ROOT / "README.md"
long_description = README.read_text(encoding="utf-8") if README.exists() else ""

INSTALL_REQUIRES = [
    "pandas==2.2.2",
    "numpy==1.26.4",
    "matplotlib",
    "igraph",
    "ete3",
    "openpyxl",
    "rich",
]

setup(
    name="daft-introgression",
    version="1.0.0",
    description=(
        "DAFT: Discordant Attachment Frequency Tool for detecting and "
        "characterizing introgression from species trees and gene trees."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Sarthak R. Mishra, Laia Pomar Pallarès, Robert Lanfear, Matthew W. Hahn",
    author_email="samishr@iu.edu",
    url="https://github.com/smishra677/DAFT",
    license="MIT",
    python_requires=">=3.9,<3.13",
    py_modules=[
        "DAFT_Test",
        "DAFT_Direction",
        "DAFT_Transform",
        "DAFT_produce_excel",
        "DAFT_produce_excel_correction",
        "clean_folder",
    ],
    packages=find_namespace_packages(
        include=[
            "daft",
            "daft.*",
            "DAFT_utils",
            "DAFT_utils.*",
        ]
    ),
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    entry_points={
        "console_scripts": [
            "daft-test=daft.cli:daft_test",
            "daft-direction=daft.cli:daft_direction",
            "daft-transform=daft.cli:daft_transform",
            "daft-excel=daft.cli:daft_excel",
            "daft-excel-correction=daft.cli:daft_excel_correction",
            "daft-clean=daft.cli:daft_clean",
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        "phylogenetics",
        "introgression",
        "gene trees",
        "hybridization",
        "discordant attachments",
        "djiNNI",
    ],
)
