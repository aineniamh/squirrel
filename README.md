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

### How it works - alignment

Squirrel maps each query genome in the input file against a reference genome specific to each clade using [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778). Using [gofasta](https://academic.oup.com/bioinformatics/article/38/16/4033/6631223), the mapping file is then converted into a multiple sequence alignment. 

For Clade II, the reference used is `NC_063383` and for Clade I, we use `NC_003310`. This means that all coordinates within an alignment will be relative to these references. A benefit of this is that within a clade, alignment files and be combined without having to recalculate the alignment. Note however that insertions relative to the reference sequence will not be included in the alignment.

Squirrel by default creates a single alignment fasta file. Using the genbank coordinates for `NC_063383` it also has the ability to extract the aligned coding sequences either as separate records or as a concatenated alignment. This can facilitate codon-aware phylogenetic or sequence analysis.

## Masking

The final alignment from squirrel has additional processing steps that includes various masking options. By default, squirrel trims the alignment to 190788 at the end of the genome to mask out one of the inverted terminal repeat (ITR) regions and pads the end of the genome with `N`. This can be disabled with the `--no-itr-mask` flag.

Squirrel performs masking (replacement with `N`) on low-complexity or repetitive regions that have been characterised for Clade I and II. These regions are defined in [to_mask.cladeii.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladeii.csv) and [to_mask.cladei.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladei.csv) respectively, and the relevant mask file will be selected depending on the `--clade` specified. This masking is on by default but can be toggled off (with `--no-mask`).

Squirrel can also accept an additional mask file in csv format (`--additional-mask`), if there are flagged sites that you wish to mask within the alignment. These sites will be masked in conjuntion with the default masking files within squirrel. The format of the mask file fits with the features format from Geneious and at a minimum should contain the following fields: "Maximum","Minimum". If squirrel is run in [QC mode](#alignment-quality-control) (`-qc` on the command like to toggle it on), it will flag sites that it believes may need to be masked from the alignment and produce a csv mask file summarising the sites. This file can then be provided to squirrel to redo the alignment with additional masking.

## Alignment Quality Control
Squirrel can run quality control (QC) on the alignment and flag certain sites to the user that may need to be masked. We recommend that the user looks at these sites in an alignment viewer to judge whether the sites should be masked or not. If `-qc` mode is toggled on, squirrel with check within the alignment for:
- <strong>Mutations that are adjacent to N bases</strong>
The rationale for this is that N sites are usually a product of low coverage regions. Mutations that occur directly adjacent to low coverage regions may be a result of mis-alignment prior to the low coverage masking and may not be real SNPs. 
- <strong>Unique mutations that clump together</strong>
If mutations are observed in only a single sequence in the genome, they are classed as unique mutations. <it>Usually</it> mutations do not clump closely together and may suggest an alignment or assembly issue. If these mutations are not shared with any other sequences, they are flagged for masking. 

## Run phylogenetic reconstruction

Signatures of APOBEC3-editing are characteristic of MPXV evolution when sustained transmission in the human population (see [O'Toole et al 2023](https://www.science.org/doi/10.1126/science.adg8116)). To also run the maximum-likelihood phylogenetics and ancestral reconstruction pipeline that characterises the APOBEC3-like and non-APOBEC3 mutations that have occurred across the evolutionary history of the virus sequences provided, run the following:

```
squirrel <your-sequences.fasta> --run-phylo --outgroups outgroup_id1,outgroup_id2
```
Note: the sequence file you provide must have the specified outgroups in it, with the IDs matching those you provide. This pipeline can accept one or more outgroup IDs. The outgroups help to root the phylogeny so that ancestral reconstruction can be performed in the appropriate direction. 

### Recommended outgroups for phylogeny mode
- <strong>Clade I</strong>
KJ642617,KJ642615,KJ642616
  
- <strong>Clade IIb</strong>
KJ642617,KJ642615

So as an example, to run phylogenetics mode on a Clade I dataset, include KJ642617, KJ642615 and KJ642616 in your test sequence file and run:

```
squirrel --clade cladei --run-phylo --outgroups KJ642617,KJ642615,KJ642616 ./test/cladeI_test.fasta
```

### How it works - phylogeny & reconstruction

Squirrel also has an optional `--run-phylo` mode that will take the newly generated alignment and build a maximum likelihood phylogeny using [IQTREE2](https://doi.org/10.1093/molbev/msaa015). It runs IQTREE ancestral state reconstruction (`-asr`), and parses the output state files, providing a branch-mapped summary of SNPs that have occurred across the phylogeny, and an output phylogeny figure with SNPs plotted along branches, coloured by whether SNPs are consistent with APOBEC3-editing or not. An outgroup (or multiple outgroups) must be specified to ensure correct rooting for the ancestral state reconstruction.

### Phylogeny-informed quality control
If you specify `-qc` when also running in phylogenetics mode (so both the `-qc` and `--run-phylo` flags), additional checks with be performed after the reconstruction. Firstly, any mutations that occur multiple times across the phylogeny (convergent mutations) are flagged for investigation, as they may highlight an issue with the tree structure or underlying sequences. 

Squirrel also flags any reversions to reference that occur in the phylogeny, which can flag issues with the assembly pipeline (often insufficient primer sequence trimming from reads or absent low-coverage masking). By default, the RefSeq records for each clade are checked against: `NC_063383` and `NC_003310`, however alternative references can be supplied with `--assembly-refs` if your sequences of interest have been assembled using a different reference, or if the primer scheme used for sequencing was constructed using a different reference.


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