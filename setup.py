import re
from setuptools import setup
  
version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('xenoGI/xenoGI.py').read(),
    re.M
    ).group(1)
 
 
with open("README.rst", "rb") as f:
    long_descr = f.read().decode("utf-8")
 
 
setup(
    name = "xenoGI",
    packages = ["xenoGI"],
    entry_points = {
        "console_scripts": ['xenoGI = xenoGI.xenoGI:main']
        },
    version = version,
    description = "Python command line application for detecting genomic island insertions in clades of microbes.",
    long_description = long_descr,
    long_description_content_type='text/x-rst',
    author = "Eliot Bush",
    author_email = "bush@hmc.edu",
    url = "https://github.com/ecbush/xenoGI",
    install_requires=['biopython','parasail'],
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    )
