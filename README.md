# squirrel

**S**ome **QUI**ck **R**earranging to **R**esolve **E**volutionary **L**inks

## Generate a quick monkeypox alignment

```
usage: squirrel <input> [options]

squirrel: Some QUIck Rearranging to Resolve Evolutionary Links

optional arguments:
  -h, --help            show this help message and exit

Input-Output options:
  input                 Input fasta file of sequences to analyse.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: current working directory
  --outfile OUTFILE     Optional output file name. Default: <input>.aln.fasta
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default: $TMPDIR
  --no-temp             Output all intermediate files, for dev purposes.

Pipeline options:
  --no-mask             Skip masking of repetitive regions. Default: masks repeat regions
  --no-itr-mask         Skip masking of end ITR. Default: masks ITR
  --extract-cds         Extract coding sequences based on coordinates in the reference
  --concatenate         Concatenate coding sequences for each genome, separated by `NNN`. Default: write out as separate records

Misc options:
  -v, --version         show program's version number and exit
  --verbose             Print lots of stuff to screen
  -t THREADS, --threads THREADS
                        Number of threads
```

## Installation


## Installation

1. Clone this repository and ``cd squirrel``
2. ``conda env create -f environment.yml``
3. ``conda activate squirrel``
4. ``pip install .``

## Check the install worked

Type (in the <strong>squirrel</strong> environment):

```
squirrel -v
```
and you should see the versions of <strong>squirrel</strong>.
