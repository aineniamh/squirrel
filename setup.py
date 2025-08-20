from setuptools import setup, find_packages
import glob
import os

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
            'tabulate',
            'baltic',
            'matplotlib',
            'mako~=1.3.10',
            'seaborn',
            'pandas~=2.3.1',
            "numpy",
            'scikit-learn==1.7.1',
            "PuLP"
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
