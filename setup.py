import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

from setuptools import setup, Extension

setup(name="ono", version="0.1",
        install_requires=[],
        scripts = ["assembly_finishing.py","assembly_finishing_objects.py"],
        ext_modules=[ \
                Extension("GSA", ["GSA.c"],  libraries=['z','divsufsort'], undef_macros=['NDEBUG'] ), \
                ],
        
        )
