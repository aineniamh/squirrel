<div align="center">

# squirrel

**S**ome **QUI**ck **R**econstruction to **R**esolve **E**volutionary **L**inks

<img src="./docs/squirrel_logo.svg" width="400">


</div>

---

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
  - [Basic Alignment](#basic-alignment)
  - [Mixed Clade Sequences](#mixed-clade-sequences)
- [Installation](#installation)
- [How It Works](#how-it-works)
  - [Alignment Process](#alignment-process)
  - [Phylogeny & Reconstruction](#phylogeny--reconstruction)
- [Masking](#masking)
  - [Default Masking Options](#default-masking-options)
---

## Quick Start

### Basic Alignment

**Align MPXV Clade II sequences:**
```bash
squirrel <your-sequences.fasta>
```

**Align MPXV Clade I sequences:**
```bash
squirrel --clade cladei <your-sequences.fasta>
```

where `<your-sequences.fasta>` is the name of your input FASTA sequence file. See [FASTA format](#fasta-format) for file format details.

### Mixed Clade Sequences

**NEW FEATURE**: If your FASTA file contains a mixture of Clade I and Clade II sequences, squirrel can now assign clade to your sequences, split them into separate alignments for Clade I and II, and run all downstream analysis per clade. A single report is produced at the end.

```bash
squirrel --clade split <your-sequences.fasta>
```

> **Note**: An EPI2ME wrapper for squirrel with the same options is available at [artic-network/squirrel-nf](https://github.com/artic-network/squirrel-nf).

---

## Installation

**Option 1: Install from bioconda (recommended)**

```bash
conda create -c bioconda -c conda-forge -n squirrel -y squirrel
conda activate squirrel
```

**Option 2: Install from source**

1. Clone this repository and `cd squirrel`
2. If `mamba` is not installed, run `conda install -n base mamba`
3. `mamba env create -f environment.yml`
4. `conda activate squirrel`
5. `pip install .`

**Verify installation:**

```bash
squirrel -v
```

You should see the version of squirrel.

---

## How It Works

### Alignment Process
**squirrel** provides a rapid way of producing reliable [alignments](#how-it-works---alignment) for MPXV and enables maximum-likelihood [phylogenetics pipeline](#phylogenetics-options-within-squirrel) tree estimation.

The reliability of tree estimation is determined by the quality of the input genome sequences. In [QC mode](#alignment-quality-control), squirrel can flag potential issues in the MPXV sequences that have been provided for alignment (e.g. SNPs near tracts of N, clusters of unique SNPs, reversions to reference alleles and convergent mutations) and outputs these in a mask file for investigation. We suggest you use this information to examine the alignment and pay close attention to the regions flagged. Squirrel can then accept this file with suggested masks and apply it to the sequences before doing phylogenetics. 

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

> *NEW FEATURE* If your FASTA file contains a mixture of Clade I and Clade II sequences, squirrel can now assign clade to your sequences, *split* them into separate alignments for Clade I and II, and run all downstream analysis per clade. A single report is produced at the end. To align a mixed FASTA file, run:
```
squirrel --clade split <your-sequences.fasta>
```

Note, an EPI2ME wrapper for squirrel with the same options is available [here](https://github.com/artic-network/squirrel-nf).

### How it works - alignment

Squirrel maps each query genome in the input file against a reference genome specific to each clade using [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778). Using [gofasta](https://academic.oup.com/bioinformatics/article/38/16/4033/6631223), the mapping file is then converted into a multiple sequence alignment. 

For Clade II, the reference used is `NC_063383` and for Clade I, we use `NC_003310`. This means that all coordinates within an alignment will be relative to these references. A benefit of this is that within a clade, alignment files and be combined without having to recalculate the alignment. Note however that insertions relative to the reference sequence will not be included in the alignment.

Squirrel by default creates a single alignment fasta file. Using the genbank coordinates for `NC_063383` it also has the ability to extract the aligned coding sequences either as separate records or as a concatenated alignment. This can facilitate codon-aware phylogenetic or sequence analysis.

If you input sequences that are entirely N or very short, squirrel may fail to align these. If all input sequences are 100% ambiguous, squirrel will exit with an error.

Some of the phylogenetic software within squirrel (iqtree, baltic) will modify or break with particular special characters in the sequence names. Squirrel will check for these special characters in the sequence names and report back if it finds any to help you curate your alignment. Special characters that we check for are: ``` ' " { ( [ ] ) } ; , ` % # : \n \t ```, however accented letters etc may also disrupt the software. For best practice, please exclude special characters from sequence names.

### Phylogeny & Reconstruction

Squirrel has the optional `--run-phylo` and `--run-apobec3-phylo` modes that will take the newly generated alignment and build a maximum likelihood phylogeny using [IQTREE2](https://doi.org/10.1093/molbev/msaa015). 

With `--run-apobec3-phylo` mode, it also runs IQTREE ancestral state reconstruction (`-asr`), and parses the output state files, providing a branch-mapped summary of SNPs that have occurred across the phylogeny, and an output phylogeny figure with SNPs plotted along branches, coloured by whether SNPs are consistent with APOBEC3-editing or not (red for APOBEC3-like mutations and yellow for other mutations). 

An outgroup (or multiple outgroups) must be specified (although this is handled internally in `--include-background` mode) to ensure correct rooting for the ancestral state reconstruction. If `--cns-qc` mode is on in conjunction with phylogenetics, reversions to reference and convergent SNPs are also flagged using the reconstruction.

---

## Masking

### Default Masking Options

The final alignment from squirrel has additional processing steps that includes various masking options:

****Sequences with a high N content**
  
  Sequences that have many ambiguous bases in them are flagged that they may want to be excluded in further analysis. This may not always be appropriate, often genomes that have a lot of ambiguity can still be informative, however if there is something unusual about a sequence, having lots of ambiguities can be a flag for wider problems (like low read count during assembly).

---

## Phylogenetics

### Maximum-Likelihood Phylogeny
**Low-complexity/Repetitive Regions** (enabled by default)
- Squirrel performs masking (replacement with `N`) on low-complexity or repetitive regions that have been characterised for Clade I and II
- These regions are defined in [to_mask.cladeii.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladeii.csv) and [to_mask.cladei.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/to_mask.cladei.csv)
- The relevant mask file will be selected depending on the `--clade` specified
- Can be disabled with `--no-mask`

### Custom Masking

**Additional Column Masking** (`--additional-mask`)
- Squirrel can accept an additional mask file in CSV format for flagged sites that you wish to mask within the alignment
- These sites will be masked in conjunction with the default masking files
- The format fits with the features format from Geneious and at a minimum should contain: "Maximum", "Minimum" fields
- If squirrel is run in [QC mode](#quality-control) (`-qc`), it will flag sites that may need to be masked and produce a CSV mask file summarising the sites

**Sequence-specific Masking** (`--sequence-mask`)
- Mask particular sites in a problematic sequence, rather than masking out the entire column across the alignment
- Requires a CSV file with a `sequence` column (must match an id in the supplied FASTA file) and a `site` column (1-based position)
- If a sequence needs multiple sites masked, specify them on separate lines with the same sequence name

---

## Quality Control
Squirrel can run quality control (QC) on the alignment and flag certain sites to the user that may need to be masked. We recommend that the user looks at these sites in an alignment viewer to judge whether the sites should be masked or not. If `-qc` mode is toggled on, squirrel with check within the alignment for:
- <strong>Mutations that are adjacent to N bases</strong>
The rationale for this is that N sites are usually a product of low coverage regions. Mutations that occur directly adjacent to low coverage regions may be a result of mis-alignment prior to the low coverage masking and may not be real SNPs. In squirrel, non-majority alleles that are present next to an N are flagged as potential sites for masking. 
- <strong>Unique mutations that clump together</strong>
If mutations are observed in only a single sequence in the genome, they are classed as unique mutations. <it>Usually</it> mutations do not clump closely together and may suggest an alignment or assembly issue. If these mutations are not shared with any other sequences, they are flagged for masking.
**Build a phylogenetic tree with specified outgroups:**

```bash
squirrel <your-sequences.fasta> --run-phylo --outgroups outgroup_id1,outgroup_id2
```

**Build a tree with automatic background sequences (no outgroup needed):**

```bash
squirrel <your-sequences.fasta> --run-phylo --include-background
```

More details on `--include-background` can be found in the [Background Mode](#background-mode) section.

### APOBEC3 Mutation Analysis

Signatures of APOBEC3-editing are characteristic of MPXV evolution during sustained transmission in the human population (see [O'Toole et al 2023](https://www.science.org/doi/10.1126/science.adg8116)). To run the maximum-likelihood phylogenetics and ancestral reconstruction pipeline that characterises the APOBEC3-like and non-APOBEC3 mutations:

```bash
squirrel <your-sequences.fasta> --run-apobec3-phylo --outgroups outgroup_id1,outgroup_id2
```

**Requirements:**
- The sequence file you provide must have the specified outgroups in it, with the IDs matching those you provide
- This pipeline can accept one or more outgroup IDs
- The outgroups help to root the phylogeny so that ancestral reconstruction can be performed in the appropriate direction

Squirrel will produce an annotated tree file that contains information on each branch about how many APOBEC3/non-APOBEC3 mutations are present on each branch, and also on each subtree beneath that branch.

### Recommended Outgroups

**NEW FEATURE**: It is no longer necessary to manually specify an outgroup. By default the recommended outgroups below are selected and pruned off for each clade I or II specifically.

- **Clade I**: KJ642617, KJ642615
- **Clade IIb**: KJ642617, KJ642615

**Example**: To run phylogenetics mode on a Clade I dataset, include KJ642617, KJ642615 and KJ642616 in your test sequence file and run:

```bash
squirrel --clade cladei --run-phylo --outgroups KJ642617,KJ642615 ./test/cladeI_test.fasta
```

### How it works - phylogeny & reconstruction

Squirrel has the optional `--run-phylo` and `--run-apobec3-phylo` modes that will take the newly generated alignment and build a maximum likelihood phylogeny using [IQTREE2](https://doi.org/10.1093/molbev/msaa015). With  `--run-apobec3-phylo` mode, it also runs IQTREE ancestral state reconstruction (`-asr`), and parses the output state files, providing a branch-mapped summary of SNPs that have occurred across the phylogeny, and an output phylogeny figure with SNPs plotted along branches, coloured by whether SNPs are consistent with APOBEC3-editing or not (red for APOBEC3-like mutations and yellow for other mutations). An outgroup (or multiple outgroups) must be specified (although this is handled internally in `--include-background` mode) to ensure correct rooting for the ancestral state reconstruction. If `--cns-qc` mode is on in conjunction with phylogenetics, reversions to reference and convergent SNPs are also flagged using the reconstruction.

### Automatic Outgroup Selection
squirrel will automatically pull out the background sequences from the relevant clade and select the appropriate outgroup from the background set. This outgroup will be pruned from the final output tree, so will not be seen in the output file.

**Automatically selected outgroups by clade:**

| Clade | Outgroup |
|-------|----------|
| cladei | KJ642615\|human\|Nigeria\|\|1978 |
| cladeia | KJ642615\|human\|Nigeria\|\|1978 |
| cladeib | KJ642613\|human\|DRC\|Equateur\|1970-09-01 |
| cladeii | KJ642613\|human\|DRC\|Equateur\|1970-09-01 |
| cladeiia | KJ642613\|human\|DRC\|Equateur\|1970-09-01 |
| cladeiib | KJ642615\|human\|Nigeria\|\|1978 |

### Background Mode

The squirrel software has a set of publically available MPXV genome sequences that include representatives of CladeIa, CladeIb, Clade IIa and CladeIIb. The sequences, the Genbank accession numbers and their clade annotations can be found in [background_sample.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background_sample.csv) and [background.fasta](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background.fasta).

If `--include-background` is run with the squirrel command, squirrel will automatically pull out the background sequences from the relevant clade and select the appropriate outgroup from the background set. 
Squirrel includes a set of publicly available MPXV genome sequences that include representatives of Clade Ia, Clade Ib, Clade IIa and Clade IIb. The sequences, GenBank accession numbers and their clade annotations can be found in:
- [background_sample.csv](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background_sample.csv)
- [background.fasta](https://github.com/aineniamh/squirrel/blob/main/squirrel/data/background.fasta)

When `--include-background` is used, squirrel will automatically:
1. Pull out the background sequences from the relevant clade
2. Select the appropriate outgroup from the background set
3. Combine them with your input sequences
4. Prune the outgroup from the final output tree

**Example**: Running `squirrel --clade cladei --include-background <sequences.fasta>` will combine all Clade I background sequences with your input file and include a Clade II outgroup sequence to correctly root the tree.

The outgroups automatically selected are the same as listed in the [Automatic Outgroup Selection](#automatic-outgroup-selection) table above.

> **Note**: If outgroups are specified in the supplied input fboth the `-qc` and `--run-phylo` flags), additional checks will be performed after the reconstruction:

**Convergent Mutations**
- Any mutations that occur multiple times across the phylogeny are flagged for investigation
- May highlight an issue with the tree structure or underlying sequences

**Reversions to Reference**
There are several options for customizing the final APOBEC3-reconstruction tree figure output:

**Re-render Tree Figure Only** (`-tfig/--tree-figure-only`)
- Rather than rerunning the entire analysis, you can re-render just the tree figure
- Requires:
  - Tree file from previous squirrel run (`-tf/--tree-file`)
  - Branch reconstruction file (`-brf/--branch-reconstruction-file`)

**Size Customization**
- `--fig-height` - Specify custom height for the final tree visualization
- `--fig-width` - Specify custom width for the final tree visualization

**Point Style Options**
- `--point-style` - Choose how mutations are represented:
**NEW FEATURE**: Squirrel can now accept a precomputed treefile, alignment and state file generated by IQTREE2 to run the APOBEC3 reconstruction.

**Note**: If the tree file has multiple likelihood metrics, it cannot be parsed by baltic.

**Required flags:**

| Flag | Description |
|------|-------------|
| `-at, --asr-tree` | Precomputed tree with ancestral states reconstructed |
| `-as, --asr-state` | Precomputed ancestral states file |
| `-aln, --asr-alignment` | Precomputed alignment file |

---

## Output Files

### Core Output Files

- **`sequences.aln.fasta`** - The alignment file, scaffolded against a clade-specific reference. By default one of the ITR regions and a curated set of problematic regions is masked as Ns.

- **`sequences.aln.tree`** - The output maximum likelihood tree file from IQTREE2 with Node labels that correspond to the reconstruction Node labels. Can be viewed in tree viewers like FigTree.

- **`sequences.aln.annotated.tree`** - The output maximum likelihood tree file from IQTREE2 with annotations on branches describing the APOBEC3-reconstruction results. Can be viewed in tree viewers like FigTree.

- **`sequences.aln.report.html`** - Summary report of analysis run.

### Reconstruction Files

- **`sequences.aln.tree.state`** and **`sequences.aln.tree.state_differences.csv`** - The output ancestral state reconstruction file from IQTREE2 and the compiled list of unambiguously variable sites from squirrel.

- **`sequences.aln.tree.branch_snps.reconstruction.csv`** - A report of individual site changes mapped to specific branches and their dinucleotide context.

- **`sequences.aln.tree.amino_acid.reconstruction.csv`** - A report of each mutation that occurs across the phylogeny: location, dinucleotide context, APOBEC3 status, gene, codon position, amino acid change, and Grantham score.

### Visualisation Files

- **`sequences.aln.tree.png`** and **`sequences.aln.tree.svg`** - Visualization of reconstructed tree showing whether mutations are consistent with APOBEC3 editing or not (red for APOBEC3-like mutations, yellow for other mutations).

- **`*tree.html`** - Interactive tree file (when using `--interactive-tree` flag).

---


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
  -o, --outdir OUTDIR   Output directory. Default: current working directory
  --outfile OUTFILE     Optional output file name. Default: <input>.aln.fasta
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default: $TMPDIR
  --no-temp             Output all intermediate files, for dev purposes.
  --epi2me-outdir EPI2ME_OUTDIR
                        Output directory for report paths in epi2me.

Alignment options:
  -qc, --seq-qc         Flag potentially problematic SNPs and sequences. Default: don't run QC
  --assembly-refs ASSEMBLY_REFS
                        References to check for `calls to reference` against.
  --no-mask             Skip masking of repetitive regions. Default: masks repeat regions
  --no-itr-mask         Skip masking of end ITR. Default: masks ITR
  --additional-mask ADDITIONAL_MASK
                        Masking additional sites provided as a csv. Needs columns `Maximum` and
                        `Minimum` in 1-base.
  --sequence-mask SEQUENCE_MASK
                        Mask sites in specific sequences in the alignment as a csv, rather than
                        the whole alignment column. Needs `sequence` and `site` (1-based)
                        column.
  -ex, --exclude EXCLUDE
                        Supply a csv file listing sequences that should be excluded from the
                        analysis.
  --extract-cds         Extract coding sequences based on coordinates in the reference
  --concatenate         Concatenate coding sequences for each genome, separated by `NNN`.
                        Default: write out as separate records
  --clade CLADE         Specify whether the alignment is primarily for `cladei` or `cladeii`
                        (can also specify a or b, e.g. `cladeia`, `cladeiib`). `split` will
                        separate input sequences into separate cladei and cladeii alignments.
                        This will determine reference used for alignment, mask file and
                        background set used if `--include-background` flag used in conjunction
                        with the `--run-phylo` option. Default: `cladeii`

Phylo options:
  -p, --run-phylo       Run phylogenetics pipeline
  -a, --run-apobec3-phylo
                        Run phylogenetics & APOBEC3-mutation reconstruction pipeline
  --outgroups OUTGROUPS
                        Specify which MPXV outgroup(s) in the alignment to use in the phylogeny.
                        These will get pruned out from the final tree.
  -bg, --include-background
                        Include a default background set of sequences for the phylogenetics
                        pipeline. The set will be determined by the `--clade` specified.
  -bf, --background-file BACKGROUND_FILE
                        Include this additional FASTA file as background to the phylogenetics.
  -bm, --binary-partition-mask
                        Calculate and write binary partition mask
  --bm-separate-dimers  Write partition mask with 0 for non-apo, 1 for GA and 2 for TC target
                        sites
  -at, --asr-tree ASR_TREE
                        Precomputed tree with ancestral states reconstructed
  -as, --asr-state ASR_STATE
                        Precomputed ancestral states file
  -aln, --asr-alignment ASR_ALIGNMENT
                        Precomputed alignment file

Tree figure options:
  -tfig, --tree-figure-only
                        Re-render tree figure custom height and width arguments. Requires: tree
                        file, branch reconstruction file, height, width.
  -tf, --tree-file TREE_FILE
                        Tree for re-rendering the figure.
  -brf, --branch-reconstruction-file BRANCH_RECONSTRUCTION_FILE
                        Reconstruction file for re-rendering the figure.
  --fig-height FIG_HEIGHT
                        Overwrite tree figure default height.
  --fig-width FIG_WIDTH
                        Overwrite tree figure default width.
  --point-style POINT_STYLE
                        Shape of points for apobec3 reconstruction figure. Options: circle,
                        square. Default: circle
  --point-justify POINT_JUSTIFY
                        Justification of points for apobec3 reconstruction figure. Options:
                        left, right. Default: left
  --interactive-tree    Create a separate interactive file built from R functions

Misc options:
  -v, --version         show program's version number and exit
  --verbose             Print lots of stuff to screen
  -t, --threads THREADS
                        Number of threads
```

### Output files

- `sequences.aln.fasta`
  
The alignment file, with alignment scaffolded against a clade-specific reference. By default one of the ITR regions and a curated set of problematic regions is masked as Ns.

- `sequences.aln.tree`

The output maximum likelihood tree file from IQTREE2 with Node labels that correspond to the reconstruction Node labels. This tree can be viewed in various tree viewers, for example FigTree. 

- `sequences.aln.annotated.tree`

The output maximum likelihood tree file from IQTREE2 with annotations on branches that describe the APOBEC3-reconstrucion results. This tree can be viewed in various tree viewers, for example FigTree. 


- `sequences.aln.tree.state` and `sequences.aln.tree.state_differences.csv`

The output ancestral state reconstruction file from IQTREE2 and the compiled list of unambiguously variable sites from squirrel.

- `sequences.aln.tree.branch_snps.reconstruction.csv`

A report of individual site changes mapped to specific branches and their dinucleotide context.

- `sequences.aln.tree.amino_acid.reconstruction.csv`

A report of each mutation that occurs across the phylogeny, their location, dinucleotide context, APOBEC3 status, which gene they're present in, codon position, amino acid change and a prediction of how extreme that amino acid change is with Grantham score.

- `sequences.aln.tree.png` and `sequences.aln.tree.svg`

Visualisation of reconstructed tree showing whether mutations are consistent with APOBEC3 editing or not.

- `sequences.aln.report.html`

Summary report of analysis run.

### What is a FASTA file <a name="fasta"></a>
<a name="fasta"></a>

A FASTA-formatted file contains sequence records. A record minimally contains two pieces of information:
1. The sequence ID (e.g. `sequence1`)
2. The sequence itself (e.g. `CGATCGAT...ACTGACT`)

**Format details:**
- The sequence ID is stored in the header line, denoted by a `>` symbol
- The header line may also contain additional information (sequence description) after the first space
- **Important**: The sequence ID must not contain whitespace (spaces or tabs)
- The sequence is stored on the following line(s)
- Sequences can be split across multiple lines for readability
- The next record does not start until the next line that begins with `>`

**Example FASTA with three sequence records:**
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
