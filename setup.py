#!/usr/bin/env python3

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))

setup(name='hivclustering',
      version="1.2.0",
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
    ],
    dependency_links = ['git+git://github.com/veg/hyphy-python.git@0.1.1#egg=HyPhy-0.1.1',
                        'git+git://github.com/veg/BioExt.git@0.17.2#egg=BioExt-0.17.3',
                        'git+git://github.com/veg/hppy.git@0.9.6#egg=hppy-0.9.6'
                       ],
    install_requires=[
        'BioExt >= 0.17.2',
        'HyPhy >= 0.1.1',
        'hppy >= 0.9.6',
        ],
     )
