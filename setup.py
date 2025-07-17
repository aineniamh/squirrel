from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from squirrel import __version__, _program

setup(name='squirrel',
      version=__version__,
      packages=find_packages(),
      scripts=[
            'squirrel/scripts/msa.smk',
            'squirrel/scripts/phylo.smk',
            'squirrel/scripts/clade.smk',
            'squirrel/scripts/reconstruction.smk',
            'squirrel/scripts/interactive_tree.R'
      ],
      package_data={"squirrel":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            'tabulate==0.8.10',
            'baltic',
            'matplotlib>=3.3.1',
            'mako',
            'seaborn',
            'pandas>2.2.0',
            "numpy>2.0.0",
            'scikit-learn==1.2.2',
            "PuLP>=2"
        ],
      description='Some QUIck Reconstruction to Resolve Evolutionary Links',
      url='https://github.com/aineniamh/squirrel',
      author='Aine OToole',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = squirrel.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
