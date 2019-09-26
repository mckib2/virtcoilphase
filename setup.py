'''Setup.py'''

from distutils.core import setup
from setuptools import find_packages

setup(
    name='virtcoilphase',
    version='0.1.0',
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/virtcoilphase',
    license='GPLv3',
    description=(
        'Virtual reference coil method to determine phase '
        'distribution'),
    long_description=open('README.rst').read(),
    install_requires=[
        "numpy>=1.17.2",
        "scipy>=1.3.1",
        "matplotlib>=3.1.1",
        "tqdm>=4.36.1"
    ],
    python_requires='>=3.6'
)
