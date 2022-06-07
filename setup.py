from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from alignHPXV import __version__, _program

setup(name='alignHPXV',
      version=__version__,
      packages=find_packages(),
      scripts=[
            'alignHPXV/scripts/msa.smk'
                ],
      package_data={"alignHPXV":["data/*"]},
      install_requires=[
            "biopython>=1.70"
        ],
      description='',
      url='https://github.com/cov-lineages/alignHPXV',
      author='Aine OToole',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = alignHPXV.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
