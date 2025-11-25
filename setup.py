from setuptools import setup, find_packages

setup(name='squirrel',
      packages=find_packages(),
      scripts=[
            'squirrel/scripts/msa.smk',
            'squirrel/scripts/phylo.smk',
            'squirrel/scripts/clade.smk',
            'squirrel/scripts/reconstruction.smk',
            'squirrel/scripts/interactive_tree.R'
      ],
      package_data={"squirrel":["data/*"]}
      )
