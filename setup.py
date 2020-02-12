#!/usr/bin/env python3

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))

setup(name='hivclustering',
      version="1.5.2",
      description='HIV molecular clustering tools',
      author='Sergei Kosakovsky Pond',
      author_email='spond@ucsd.edu',
      url='http://github.com/veg/hivclustering',
      license='MIT License',
      packages=['hivclustering'],
      package_data={'hivclustering': [
            'data/HBL/*.bf',
    ]},
    scripts=[
        'scripts/hivnetworkcsv',
        'scripts/TNS',
        'scripts/hivnetworkannotate'
    ],
    install_requires=[
        'bioext >= 0.19.0',
        'hyphy-python >= 0.1.9',
        'hppy >= 0.9.9',
        ],
     )
