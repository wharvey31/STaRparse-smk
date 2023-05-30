import json
import pandas as pd
import glob
import os, sys
import pyvcf


configfile: "config.yaml"


REFERENCE = (config["REF"],)
VARIANT_CATALOG = (config["VARIANT_CATALOG"],)
MANIFEST = config.get("MANIFEST", "manifest.tab")
GENOME_BUILD = config.get("GENOME_BUILD", "GRCh38")

SDIR = os.path.realpath(os.path.dirname(srcdir("Snakefile")))

manifest_df = pd.read_csv(MANIFEST, sep="\t")


def find_aln(wildcards):
    return manifest_df.at[wildcards.sample, "ALN"]


def find_sex(wildcards):
    return manifest_df.at[wildcards.sample, "SEX"]


wildcard_constraints:
    sample="|".join(manifest_df.index),


localrules:
    all,


rule all:
    input:
        expand(
            "results/vcf_parse/CGG_Repeats_{sample}_filtered.tsv",
            sample=manifest_df.index,
        ),


rule expansion_hunter:
    input:
        aln=find_aln,
        catalog=VARIANT_CATALOG,
        ref=REFERENCE,
    output:
        json="results/eh/{sample}.json",
        vcf="results/eh/{sample}.vcf",
    params:
        sex=find_sex,
    conda:
        "envs/eh.yaml"
    threads: 1
    resources:
        mem=8,
        hrs=24,
    shell:
        """
        ExpansionHunter --reads {input.aln} --reference {input.ref} --varaint-catalog {input.catalog} --output-prefix $( echo {input.vcf} | sed 's/.vcf//' ) --sex {params.sex}
        """


rule vcfparse:
    input:
        vcf=rules.expansion_hunter.output.vcf,
    output:
        csv="tmp/vcf_parse/CGG_Repeats_{sample}.tsv",
    threads: 1
    resources:
        mem=8,
        hrs=24,
    run:
        with open(output.csv, "w") as outfile:
            header = "Call,Sample_ID,Chr,Start,End,GT,Ref_Units,Allele1_Units,Allele2_Units".replace(
                ",", "\t"
            )
            outfile.write(f"{header}\n")
            my_vcf = vcf.Reader(filename=input.vcf)
            for record in my_vcf:
                samp = record.samples
                for i in samp:
                    if i["GT"] == ".":
                        continue
                    Call = record.INFO["REPID"]
                    chrom = record.CHROM
                    start = int(record.POS)
                    end = int(record.INFO["END"])
                    units = i["REPCN"]
                    gt = i["GT"]
                    ref = round(record.INFO["REF"])
                    if "/" not in units:
                        Allele1_Units = units
                        Allele2_Units = 0
                    else:
                        Allele1_Units = units.split("/")[0]
                        Allele2_Units = units.split("/")[1]
                    strings = [
                        str(val)
                        for val in [
                            Call,
                            base,
                            chrom,
                            start,
                            end,
                            gt,
                            ref,
                            Allele1_Units,
                            Allele2_Units,
                        ]
                    ]
                    strings = "\t".join(strings)
                    outfile.write(f"{strings}\n")


rule jsonparse:
    input:
        json=rules.expansion_hunter.output.json,
    output:
        csv="tmp/json_parse/ReadCoverage_{sample}.tsv",
    threads: 1
    resources:
        mem=8,
        hrs=24,
    run:
        #    COLLECT JSON FILES FROM SPECIFIED INPUT PATH
        #    INITIALIZE STOREAGE FILE
        with open(output.csv, "w+") as outfile:
            header = "Sample_ID,Call_ID,MaxSpanningRead,MaxFlankingRead,MaxInrepeatRead".replace(
                ",", "\t"
            )
            outfile.write(f"{header}\n")
            #    LOOP THROUGH ALL JSON FILES AT GIVEN DIRECTORY
            with open(j) as json_file:
                data = json.load(json_file)["LocusResults"]
            #    LOOP THROUGH EACH READ WITHIN THE CURRENT JSON FILE
            for x in data:
                repeat = data[x]
                #    ENSURE THAT READ WAS DECTECTED AND DATA IS PRESENT
                if len(repeat) < 5:
                    continue
                else:
                    #    COLLECT THE LENGTH AND NUMBER OF FLANKING, SPANNING, AND INREPEAT READS
                    span = str(repeat["Variants"][x]["CountsOfSpanningReads"])
                    flank = str(repeat["Variants"][x]["CountsOfFlankingReads"])
                    inread = str(repeat["Variants"][x]["CountsOfInrepeatReads"])

                    strings = [base, x, span, flank, inread]
                    strings = "\t".join(strings)
                    outfile.write(f"{strings}\n")


rule filter_r:
    input:
        reads=rules.vcfparse.output.csv,
        coverage=rules.jsonparse.output.csv,
    output:
        csv="results/vcf_parse/CGG_Repeats_{sample}_filtered.tsv",
    params:
        script_dir=os.path.join(SDIR, "scripts/"),
    conda:
        "envs/r_conda.yaml"
    threads: 1
    resources:
        mem=8,
        hrs=24,
    script:
        """
        scripts/filter.R
        """


rule summaries:
    input:
        csv=rules.filter_r.output.csv,
    output:
        csv="results/filter/{sample}/{sample}_by_locus.tsv",
    params:
        build=GENOME_BUILD,
        script_dir=os.path.join(SDIR, "scripts/"),
    conda:
        "envs/r_conda.yaml"
    threads: 1
    resources:
        mem=8,
        hrs=24,
    script:
        """
        scripts/summaries.R
        """
