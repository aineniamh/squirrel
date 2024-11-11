# squirrel

**S**ome **QUI**ck **R**econstruction to **R**esolve **E**volutionary **L**inks

Why use squirrel? 

The MPXV genome is pretty challenging to work with and do reliable phylogenetics on. It is large (~200kb), has tracts of low complexity and repetitive regions, and has large deletions, which can lead to difficulties producing a reliable alignment. With squirrel, we provide a rapid way of producing reliable [alignments](#how-it-works---alignment) for MPXV and also enable maximum-likelihood [phylogenetics pipeline](#phylogenetics-options-within-squirrel) tree estimation. 

The reliability of tree estimation is determined by the quality of the input genome sequences. In [QC mode](#alignment-quality-control), squirrel can flag potential issues in the MPXV sequences that have been provided for alignment (e.g. SNPS near tracts of N, clusters of unique SNPs, reversions to reference alleles and convergent mutations) and outputs these in a mask file for investigation. We suggest you use this information to examine the alignment and pay close attention to the regions flagged. Squirrel can then accept this file with suggested masks and apply it to the sequences before doing phylogenetics. 

Enrichment of APOBEC3-mutations in the MPXV population are a signature of sustained human-to-human transmission. Identifying APOBEC3-like mutations in MPXV genomes from samples in a new outbreak can be a piece of evidence to support sustained human transmission of mpox. Squirrel can run an [APOBEC3-reconstruction](#phylogenetics-options-within-squirrel) and map these mutations onto the phylogeny.

Squirrel produces a HTML report summarising some of the key outputs.

## Generate a quick MPXV alignment
> To align MPXV Clade II sequences, run:
```
squirrel <your-sequences.fasta>
```

>To align MPXV Clade I sequences, run:

```
squirrel --clade cladei <your-sequences.fasta>
```
where `<your-sequences.fasta>` is the name of your input FASTA sequence file. Click [here](#fasta) see what a FASTA formatted file looks like.

Note, an EPI2ME wrapper for squirrel with the same options is available [here](https://github.com/artic-network/squirrel-nf).

### How it works - alignment

Squirrel maps each query genome in the input file against a reference genome specific to each clade using [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778). Using [gofasta](https://academic.oup.com/bioinformatics/article/38/16/4033/6631223), the mapping file is then converted into a multiple sequence alignment. 

For Clade II, the reference used is `NC_063383` and for Clade I, we use `NC_003310`. This means that all coordinates within an alignment will be relative to these references. A benefit of this is that within a clade, alignment files and be combined without having to recalculate the alignment. Note however that insertions relative to the reference sequence will not be included in the alignment.

Squirrel by default creates a single alignment fasta file. Using the genbank coordinates for `NC_063383` it also has the ability to extract the aligned coding sequences either as separate records or as a concatenated alignment. This can facilitate codon-aware phylogenetic or sequence analysis.

## Masking

The final alignment from squirrel has additional processing steps that includes various masking options. By default, squirrel trims the alignment to 190788 at the end of the genome to mask out one of the inverted terminal repeat (ITR) regions and pads the end of the genome with `N`. This can be disabled with the `--no-itr-mask` flag.

Squirrel performs masking (replacement with `N`) on low-complexity or repetitive regions that have been characterised for Clade I and II. These regions are defined in [to_mask.cladeii.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladeii.csv) and [to_mask.cladei.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladei.csv) respectively, and the relevant mask file will be selected depending on the `--clade` specified. This masking is on by default but can be toggled off (with `--no-mask`).

Squirrel can also accept an additional mask file in csv format (`--additional-mask`), if there are flagged sites that you wish to mask within the alignment. These sites will be masked in conjuntion with the default masking files within squirrel. The format of the mask file fits with the features format from Geneious and at a minimum should contain the following fields: "Maximum","Minimum". If squirrel is run in [QC mode](#alignment-quality-control) (`-qc` on the command like to toggle it on), it will flag sites that it believes may need to be masked from the alignment and produce a csv mask file summarising the sites. This file can then be provided to squirrel to redo the alignment with additional masking.

Squirrel now offers a sequence-specific mask option for making particular sites in a problematic sequence, rather than masking out the entire column across the alignment. This can be achieved with `--sequence-mask` and requires the user to supply a csv file with a `sequence` column that contains the sequence name (must match an id in the supplied FASTA file) and a `site` column, with a 1-based position to mask. If a sequence needs multiple sites masked, specify them on separate lines in the file with the same sequence name.

## Alignment Quality Control
Squirrel can run quality control (QC) on the alignment and flag certain sites to the user that may need to be masked. We recommend that the user looks at these sites in an alignment viewer to judge whether the sites should be masked or not. If `-qc` mode is toggled on, squirrel with check within the alignment for:
- <strong>Mutations that are adjacent to N bases</strong>
The rationale for this is that N sites are usually a product of low coverage regions. Mutations that occur directly adjacent to low coverage regions may be a result of mis-alignment prior to the low coverage masking and may not be real SNPs. In squirrel, non-majority alleles that are present next to an N are flagged as potential sites for masking. 
- <strong>Unique mutations that clump together</strong>
If mutations are observed in only a single sequence in the genome, they are classed as unique mutations. <it>Usually</it> mutations do not clump closely together and may suggest an alignment or assembly issue. If these mutations are not shared with any other sequences, they are flagged for masking.
- <strong>Sequences with a high N content</strong>
Sequences that have many ambiguous bases in them are flagged that they may want to be excluded in further analysis. This may not always be appropriate, often genomes that have a lot of ambiguity can still be informative, however if there is something unusual about a sequence, having lots of ambiguities can be a flag for wider problems (like low read count during assembly).

## Phylogenetics options within squirrel

To build a maximum-likelihood phylogeny with [IQTREE2](https://doi.org/10.1093/molbev/msaa015) with the alignment generated, run the following:

```
squirrel <your-sequences.fasta> --run-phylo --outgroups outgroup_id1,outgroup_id2
```

If you wish to include a background set of sequences and not have to specify an outgroup, run the following:

```
squirrel <your-sequences.fasta> --run-phylo --include-background
```

More details on `--include-background` can be found [here](#automatically-include-background-and-outgroups-with-include-background-mode).



Signatures of APOBEC3-editing are characteristic of MPXV evolution when sustained transmission in the human population (see [O'Toole et al 2023](https://www.science.org/doi/10.1126/science.adg8116)). To also run the maximum-likelihood phylogenetics and ancestral reconstruction pipeline that characterises the APOBEC3-like and non-APOBEC3 mutations that have occurred across the evolutionary history of the virus sequences provided, run the following:

```
squirrel <your-sequences.fasta> --run-apobec3-phylo --outgroups outgroup_id1,outgroup_id2
```
Note: the sequence file you provide must have the specified outgroups in it, with the IDs matching those you provide. This pipeline can accept one or more outgroup IDs. The outgroups help to root the phylogeny so that ancestral reconstruction can be performed in the appropriate direction. 

### Recommended outgroups for phylogeny mode
- <strong>Clade I</strong>
KJ642617,KJ642615
  
- <strong>Clade IIb</strong>
KJ642617,KJ642615

So as an example, to run phylogenetics mode on a Clade I dataset, include KJ642617, KJ642615 and KJ642616 in your test sequence file and run:

```
squirrel --clade cladei --run-phylo --outgroups KJ642617,KJ642615 ./test/cladeI_test.fasta
```


### How it works - phylogeny & reconstruction

Squirrel has the optional `--run-phylo` and `--run-apobec3-phylo` modes that will take the newly generated alignment and build a maximum likelihood phylogeny using [IQTREE2](https://doi.org/10.1093/molbev/msaa015). With  `--run-apobec3-phylo` mode, it also runs IQTREE ancestral state reconstruction (`-asr`), and parses the output state files, providing a branch-mapped summary of SNPs that have occurred across the phylogeny, and an output phylogeny figure with SNPs plotted along branches, coloured by whether SNPs are consistent with APOBEC3-editing or not. An outgroup (or multiple outgroups) must be specified (although this is handled internally in `--include-background` mode) to ensure correct rooting for the ancestral state reconstruction. If `--cns-qc` mode is on in conjunction with phylogenetics, reversions to reference and convergent SNPs are also flagged using the reconstruction.

### Automatically include background and outgroups with include-background mode 

The squirrel software has a set of publically available MPXV genome sequences that include representatives of CladeIa, CladeIb, Clade IIa and CladeIIb. The sequences, the Genbank accession numbers and their clade annotations can be found in [background_sample.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background_sample.csv) and [background.fasta](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background.fasta).

If `--include-background` is run with the squirrel command, squirrel will automatically pull out the background sequences from the relevant clade and select the appropriate outgroup from the background set. 

For example,  if `squirrel --clade cladei --include-background <sequences.fasta>` is run, squirrel will combine all Clade I background sequences with the input `<sequences.fasta>` file and also include a Clade II outgroup sequence to correctly root the tree. This outgroup will be pruned from the final output tree, so will not be seen in the output file.

The respective outgroups automatically selected from the background set are:
```
    cladei: "KJ642615|human|Nigeria||1978"
    cladeia: "KJ642615|human|Nigeria||1978"
    cladeib: "KJ642613|human|DRC|Equateur|1970-09-01"
    cladeii: "KJ642613|human|DRC|Equateur|1970-09-01"
    cladeiia: "KJ642613|human|DRC|Equateur|1970-09-01"
    cladeiib: "KJ642615|human|Nigeria||1978"
```
>Note that in this mode, if outgroups are specified in the supplied input file, they will not be used and they will remain in the final tree.

### Phylogeny-informed quality control
If you specify `-qc` when also running in phylogenetics mode (so both the `-qc` and `--run-phylo` flags), additional checks with be performed after the reconstruction. Firstly, any mutations that occur multiple times across the phylogeny (convergent mutations) are flagged for investigation, as they may highlight an issue with the tree structure or underlying sequences. 

Squirrel also flags any reversions to reference that occur in the phylogeny, which can flag issues with the assembly pipeline (often insufficient primer sequence trimming from reads or absent low-coverage masking). By default, the RefSeq records for each clade are checked against: `NC_063383` and `NC_003310`, however alternative references can be supplied with `--assembly-refs` if your sequences of interest have been assembled using a different reference, or if the primer scheme used for sequencing was constructed using a different reference.

### APOBEC3-reconstruction tree figure customisation

Occasionally, you may want to adjust the final APOBEC3-reconstruction tree figure output. There are now a number of options to help you achieve this. 

Firstly, rather than having to rerun the entire analysis, there is now an option that allows the user to re-render the tree figure only (`-tfig/--tree-figure-only`). The user must supply a tree file and matching branch reconstruction file that has been output by a previous squirrel run. This can be done with the `-brf/--branch-reconstruction-file` and `-tf/--tree-file` arguments. 

Alongside these files, the user can then specify a custom height (`--fig-height`) or width (`--fig-width`) for the final tree vizualisation.

The user may now also specify whether the reconstructed mutations vizualised on the branch are either represented by a `circle` hovering over the branch or a `square` spanning the branch with the `--point-style` argument, and whether they want the points to begin stacking from the `left` or `right` with `--point-justify`. 

## Installation

Install from bioconda with `conda` or `mamba`, e.g. `conda create -c bioconda -c conda-forge -n squirrel -y squirrel`.

Or:

1. Clone this repository and ``cd squirrel``
2. if `mamba` is not installed, run `conda install -n base mamba`
3. ``mamba env create -f environment.yml``
4. ``conda activate squirrel``
5. ``pip install .``

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

options:
  -h, --help            show this help message and exit

Input-Output options:
  input                 Input fasta file of sequences to analyse.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: current working directory
  --outfile OUTFILE     Optional output file name. Default: <input>.aln.fasta
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default: $TMPDIR
  --no-temp             Output all intermediate files, for dev purposes.

Pipeline options:
  -qc, --seq-qc         Flag potentially problematic SNPs and sequences. Note that this will also run phylo mode, so you will need to specify both outgroup sequences and provide an assembly
                        reference file. Default: don't run QC
  --assembly-refs ASSEMBLY_REFS
                        References to check for `calls to reference` against.
  --no-mask             Skip masking of repetitive regions. Default: masks repeat regions
  --no-itr-mask         Skip masking of end ITR. Default: masks ITR
  --additional-mask ADDITIONAL_MASK
                        Masking additional sites provided.
  --extract-cds         Extract coding sequences based on coordinates in the reference
  --concatenate         Concatenate coding sequences for each genome, separated by `NNN`. Default: write out as separate records
  --clade CLADE         Specify whether the alignment is primarily for `cladei` or `cladeii` (can also specify a or b, e.g. `cladeia`, `cladeiib`). This will determine reference used for
                        alignment, mask file and background set used if `--include-background` flag used in conjunction with the `--run-phylo` option. Default: `cladeii`
  -p, --run-phylo       Run phylogenetics pipeline
  -a, --run-apobec3-phylo
                        Run phylogenetics & APOBEC3-mutation reconstruction pipeline
  --outgroups OUTGROUPS
                        Specify which MPXV outgroup(s) in the alignment to use in the phylogeny. These will get pruned out from the final tree.
  -bg, --include-background
                        Include a default background set of sequences for the phylogenetics pipeline. The set will be determined by the `--clade` specified.
  -bf BACKGROUND_FILE, --background-file BACKGROUND_FILE
                        Include this additional FASTA file as background to the phylogenetics.

Misc options:
  -v, --version         show program's version number and exit
  --verbose             Print lots of stuff to screen
  -t THREADS, --threads THREADS
                        Number of threads
```

### What is a FASTA file <a name="fasta"></a>

A FASTA-formatted file contains sequence records. A record minimally contains two pieces of information, the sequence ID (e.g. `sequence1`) and the sequence itself (e.g. `CGATCGAT...ACTGACT`). The sequence ID is stored in the header line, which is denoted by a `>` symbol. The header line may also contain additional information (called the sequence description), which can be found after the first space on the header line (e.g. `some_extra_information` in the example below). This is why it is important that the sequence ID does not contain whitespace (i.e. spaces or ` `). The sequence itself is then stored on the following line. Often the sequence is split across multiple lines for readability, but note that the next record does not start until the next line that begins with `>`. An example of this split line display is below for the record containing sequence2. 

Example FASTA with three sequence records:
```
>sequence1 some_extra_information
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
>sequence2 some_MORE_extra_information
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
GGCTAGCTAGCGTAGCTAGCGCATTACGTACTACT
TGCTAGCTAGCGTAGCTAGCGCATTACGTACTACA
ACGTAGTCATAGTCGTACTGAC
>sequence3
AGCTAGCTAGCGTAGCTAGCGCATTACGTACTACG
```
