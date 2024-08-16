# squirrel

**S**ome **QUI**ck **R**earranging to **R**esolve **E**volutionary **L**inks

## Generate a quick MPXV alignment
> To align MPXV Clade II sequences, run:
```
squirrel <your-sequences.fasta>
```

>To align MPXV Clade I sequences, run:

```
squirrel --clade cladeii <your-sequences.fasta>
```
where `<your-sequences.fasta>` is the name of your input FASTA sequence file. Click [here](#fasta) see what a FASTA formatted file looks like. 

### How it works - alignment and options



Squirrel maps each query genome in the input file against the NC_063383 reference genome using [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778). 

Using [gofasta](https://academic.oup.com/bioinformatics/article/38/16/4033/6631223), the mapping file is then converted into a multiple sequence alignment. 

It then trims to 190788 at the end of the genome to mask out one of the ITR regions and pads the end of the genome with `N`. It performs masking (replacement with `N`) on low-complexity or repetitive regions, defined in [to_mask.cladeii.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladeii.csv). The masking is on by default but can be toggled off (with `--no-mask`) and the 


Squirrel by default creates a single alignment fasta file. Using the genbank coordinates for `NC_063383` it also has the ability to extract the aligned coding sequences either as separate records or as a concatenated alignment. This can facilitate codon-aware phylogenetic or sequence analysis.

## 

## Run reconstruction

```
squirrel <your-sequences.fasta> --run-phylo --outgroups outgroup_id1,outgroup_id2
```
Note: the sequence file you provide must have the specified outgroups in it, with the IDs matching those you provide. This pipeline can accept one or more outgroup IDs.


## How it works - phylogeny & reconstruction

Squirrel also has an optional `--run-phylo` mode that will take the newly generated alignment and build a maximum likelihood phylogeny using iqtree. It runs iqtree for ancestral state reconstruction too, and parses the output files providing a branch-mapped summary of SNPs that have occurred across the phylogeny, and an output phylogeny figure with SNPs plotted along branches, coloured by whether SNPs are consistent with APOBEC3-editing or not. An outgroup (or multiple outgroups) must be specified to ensure correct rooting for the ancestral state reconstruction.

## Recommended outgroups for phylogeny mode
Clade I
- KJ642617,KJ642615,KJ642616
  
Clade IIb
- KJ642617,KJ642615

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

## Full usage


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
  -qc, --seq-qc         Flag potentially problematic SNPs and sequences. Note that this will also run phylo mode, so you will need to specify both outgroup sequences and provide an assembly reference file. Default: don't run
                        QC
  --assembly-refs ASSEMBLY_REFS
                        References to check for `calls to reference` against.
  --no-mask             Skip masking of repetitive regions. Default: masks repeat regions
  --no-itr-mask         Skip masking of end ITR. Default: masks ITR
  --additional-mask ADDITIONAL_MASK
                        Masking additional sites provided.
  --extract-cds         Extract coding sequences based on coordinates in the reference
  --concatenate         Concatenate coding sequences for each genome, separated by `NNN`. Default: write out as separate records
  --clade CLADE         Specify whether the alignment is primarily for `cladei` or `cladeii` (will determine reference used for alignment). Default: `cladeii`
  -p, --run-phylo       Run phylogenetic reconstruction pipeline
  --outgroups OUTGROUPS
                        Specify which MPXV outgroup(s) in the alignment to use in the phylogeny. These will get pruned out from the final tree.

Misc options:
  -v, --version         show program's version number and exit
  --verbose             Print lots of stuff to screen
  -t THREADS, --threads THREADS
                        Number of threads
```

### What is a FASTA file <a name="fasta"></a>

```
>sequence1 some_extra_information
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
>sequence2 some_MORE_extra_information
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
>sequence3
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
```