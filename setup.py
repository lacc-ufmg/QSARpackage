from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='qsar_package',
    version='0.0.1-dev',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'qsar_package=qsar_package.__main__:main',
        ],
    },
    install_requires=[
        "click",
        "numpy",
        # "pandas",
        # "scipy",
        # "openbabel",
        # "rdkit",
        # "openmm",
    ],
   
    author="Martins, J. P. A; Reis Filho, H. M.",
    author_email="jpam@qui.ufmg.br,helitonmrf@ufmg.br",
    description="QSAR analysis in Python",
    long_description=long_description,
    url="https://github.com/lqcapf/qsarpackage",
    # package_dir={"qsar_package": "qsar_package"},
    package_data={"": ["lqtagrid/defaultFiles/AtomProva.atp", "lqtagrid/defaultFiles/ffcargasnb.itp", "lqtagrid/defaultFiles/*.csv", "md/static/*.txt"]},
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.8',
    license_files = ('LICENSE',),
)