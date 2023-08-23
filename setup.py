#!/usr/bin/env python3

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))

setup(name='hivclustering',
      version="1.6.7",
      description='HIV molecular clustering tools',
      author='Sergei Kosakovsky Pond',
      author_email='spond@ucsd.edu',
      url='http://github.com/veg/hivclustering',
      license='MIT License',
      packages=['hivclustering'],
      package_data={'hivclustering': [
          'data/HBL/*.bf',
      ]},
      extras_require={
          'edgefiltering': ['bioext >= 0.19.0', 'hyphy-python >= 0.1.11', 'hppy >= 0.9.9'],
      },
      entry_points={
          'console_scripts': [
              'hivnetworkcsv = hivclustering.scripts:hivnetworkcsv',
              'hivnetworkannotate = hivclustering.scripts:hivnetworkannotate',
              'TNS = hivclustering.scripts:TNS',
          ]
      },
      )
