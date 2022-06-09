import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan


if config[KEY_EXTRACT_CDS]:
    rule all:
        input:
            os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILE]),
            os.path.join(config[KEY_OUTDIR],config[KEY_CDS_OUTFILE])
else:
    rule all:
        input:
            os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILE])

rule align_to_reference:
    input:
        fasta = config[KEY_INPUT_FASTA],
        reference = config[KEY_REFERENCE_FASTA]
    params:
        trim_start = 0,
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_TEMPDIR],"mapped.sam")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"msa.fasta")
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/minimap2_sam.log")
    run:
        # the first line of this streams through the fasta and replaces '-' in sequences with empty strings
        # this could be replaced by a python script later
        #  {{ gsub(" ","_",$0); }} {{ gsub(",","_",$0); }}
        shell_command =  """ | awk '{{ if ($0 !~ /^>/) {{ gsub("-", "",$0); }} print $0; }}'   | \
            awk '{{ {{ gsub(" ", "_",$0); }} {{ gsub(",", "_",$0); }} print $0; }}'  | \
            minimap2 -a -x asm20 -rmq=no --junc-bonus=0 --for-only --sam-hit-only --secondary=no --score-N=0  -t  {workflow.cores} {input.reference:q} - -o {params.sam:q} &> {log:q}
            gofasta sam toMultiAlign \
                -s {params.sam:q} \
                -t {workflow.cores} \
                --reference {input.reference:q} \
                --trimstart {params.trim_start} \
                --trimend {params.trim_end} \
                --trim \
                --pad -o '{output.fasta}' &> {log:q}
            """

        shell_command = "cat \"{input.fasta}\"" + shell_command
        shell(shell_command)

rule mask_repetitive_regions:
    input:
        mask = config[KEY_TO_MASK],
        fasta = rules.align_to_reference.output.fasta
    output:
        os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILE])
    run:
        if not config["no_mask"]:
            mask_sites = []
            total_masked = 0
            with open(input.mask, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    start = int(row["Minimum"]) - 1
                    end = int(row["Maximum"])
                    length = int(row["Length"])
                    mask_sites.append((start,end,length))
                    total_masked += length

            records = 0
            with open(output[0],"w") as fw:
                for record in SeqIO.parse(input.fasta,"fasta"):
                    records+=1
                    new_seq = str(record.seq)
                    for site in mask_sites:
                        new_seq = new_seq[:site[0]] + ("N"*site[2]) + new_seq[site[1]:]
                    n_diff = str(record.seq).count("N") - new_seq.count("N")
                    fw.write(f">{record.description}\n{new_seq}\n")
            
            print(green(f"{records} masked, aligned sequences written to: ") + f"{output[0]}")
        else:
            shell("cp {input.fasta:q} {output[0]:q}")
            print(green(f"Aligned sequences written to: ") + f"{output[0]}")

rule extract_cds:
    input:
        boundaries = config[KEY_GENE_BOUNDARIES],
        fasta = rules.mask_repetitive_regions.output[0]
    output:
        os.path.join(config[KEY_OUTDIR],config[KEY_CDS_OUTFILE])
    run:
        genes = {}
        gene_id = 0
        with open(input.boundaries, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene_id +=1
                name = f"{row['Name'].replace(' ','_')}_{gene_id}"
                start = int(row["Minimum"]) - 1
                end = int(row["Maximum"])
                length = int(row["Length"])
                direction = row["Direction"]
                # print(name,start,end,direction)
                genes[name]=(start,end,length,direction)

        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.fasta,"fasta"):
                full_genome = record.seq
                extractions = []
                for gene in genes:

                    start,end,length,direction = genes[gene]
                    extraction = full_genome[start:end]
                    if len(extraction)%3 != 0:
                        print("Extraction not a multiple of 3: ",gene,", length: ",len(extraction)," coords: ,",start,end)
                    if direction == "reverse":
                        extraction = extraction.reverse_complement()
                    extraction_name = f"{record.description}|{gene}"
                    extractions.append((extraction_name,start,end, direction, extraction))

                if config[KEY_CONCATENATE]:
                    seqs_to_write = [str(i[-1]) for i in extractions]
                    new_seq = "NNN".join(seqs_to_write)
                    fw.write(f">{record.description}\n{new_seq}\n")
                else:
                    for i in extractions:
                        fw.write(f">{i[0]}|{i[1]}-{i[2]}|{i[3]}\n{i[-1]}\n")
        print(green(f"CDS sequences written to: ") + f"{output[0]}")



