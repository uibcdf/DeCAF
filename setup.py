#!/usr/bin/env python

from setuptools import setup

setup(
    name='DeCAF',
    version='2.0.0',
    description='Discrimination, Comparison, Alignment tool for small molecules',
    author='Marta M. Stepniewska',
    author_email='martasd@ibb.waw.pl',
    url='https://bitbucket.org/marta-sd/decaf',
    license='BSD',
    packages=['decaf',
              'decaf.toolkits'],
)
