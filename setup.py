#!/usr/bin/env python

from setuptools import setup

setup(
    name='DeCAF',
    version='2.0.0',
    description='Discrimination, Comparison, Alignment tool for small molecules',
    keywords=['cheminformatics', 'small molecules', 'pharmacophore', 'graph',
              'alignment'],
    author='Marta M. Stepniewska',
    author_email='martasd@ibb.waw.pl',
    url='https://bitbucket.org/marta-sd/decaf',
    license='BSD',
    packages=['decaf',
              'decaf.toolkits',
              'tests'],
    setup_requires=['numpy>=1.8'],
    install_requires=['numpy>=1.8', 'scipy>=0.13', 'matplotlib>=1.3.1'],
    test_suite='tests',
    package_data={'tests': ['tests/*.phar']}
)
