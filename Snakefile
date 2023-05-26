import stringr
import plyr
import optparse
import reshape2
import dplyr
import plyr
import json
import pandas as pd
import glob
import os
import subprocess
import pyvcf


rule all:

rule vcfparse:
    input:
        vcf = find_vcf_input
    output:
        csv = 'tmp/vcf_parse/CGG_Repeats_{sample}.tsv'
    run:
        with open(output.csv,"w") as f:
            header = 'Call,Sample_ID,Chr,Start,End,GT,Ref_Units,Allele1_Units,Allele2_Units'.replace(',', '\t')
            f.write(f"{header}\n")
            my_vcf = vcf.Reader(filename=input.vcf)
            for record in my_vcf:
                samp = record.samples
                for i in samp:
                     if i['GT'] == ".":
                         continue
                     Call = record.INFO['REPID']
                     chrom = record.CHROM
                     start = int(record.POS)
                     end = int(record.INFO["END"])
                     units = i['REPCN']
                     gt = i['GT']
                     ref = round(record.INFO["REF"])
                     if '/' not in units:
                         Allele1_Units=units
                         Allele2_Units=0
                     else:
                         Allele1_Units=units.split("/")[0]
                         Allele2_Units=units.split("/")[1]
                     strings = [str(val) for val in [Call , base, chrom, start, end, gt, ref, Allele1_Units, Allele2_Units] ]
                     strings="\t".join(strings)
                     f.write(f"{strings}\n")

rule jsonparse:
    input:
        json = find_json_input
    output:
        csv = 'tmp/json_parse/ReadCoverage_{sample}.tsv'
    run:
        #    COLLECT JSON FILES FROM SPECIFIED INPUT PATH
        #    INITIALIZE STOREAGE FILE
        with open(output.csv, 'w+') as f:\
            header = "Sample_ID,Call_ID,MaxSpanningRead,MaxFlankingRead,MaxInrepeatRead".replace(',', '\t')
            f.write(f"{header}\n")
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
                    span = str(repeat['Variants'][x]['CountsOfSpanningReads'])
                    flank = str(repeat['Variants'][x]['CountsOfFlankingReads'])
                    inread = str(repeat['Variants'][x]['CountsOfInrepeatReads'])

                    strings = [base, x, span, flank, inread]
                    strings = '\t'.join(strings)
                    f.write(f"{strings}\n")

rule filter_r:
    input:
        reads = rules.vcfparse.output.csv,
        coverage = rules.jsonparse.output.csv
    output:
        csv = 'tmp/vcf_parse/CGG_Repeats_{sample}_filtered.tsv',
    params:
        script_dir = os.path.join(SDIR, "scripts/")
    script:
        """
        scripts/filter.R
        """

rule summaries:
    input:
        csv = rules.filter_r.output.csv,
    output:
        csv = "temp/{sample}_by_locus.csv"
    params:
        build = GENOME_BUILD,
        script_dir = os.path.join(SDIR, "scripts/")
    script:
        """
        scripts/script.R
        """




