# tWPS (tx) below were adapted from scripts in cfDNA (2016, Cell, Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin) 
# 174:         awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ print $1,$2-{params.flank},$3+{params.flank},$4,$5,$6 }}' |
# test add shell.prefix due to plot_overlay.py error
#shell.prefix("set -x; set -e;")

from snakemake.utils import validate
import pandas as pd
from os.path import exists

#configfile: "config/config.yml"

validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="schemas/regions.schema.yaml")


def get_WPS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_WPS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_COV.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_STARTS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_WPS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_WPS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_COV.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/background/{{GENOME}}/table/{{target_region}}--{ref_SAMPLE}_STARTS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )

#only use first row lenth, treat all other regions has the same length
def get_length(input):
    if exists(input):
        df = pd.read_csv(input, sep="\t", header=None)
        length = df[2] - df[1]
        return length[0]
    else:
        length=2000
        #print(f"{input} does not exist: using length of {length}")
        return length
	# length: may used for enough calculating extended region ;only return 1st value ?


rule all:
    input:
        # expand(
        #     expand(
        #         "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        #         zip,
        #         ID=samples["ID"],
        #         GENOME=samples["genome_build"],
        #         allow_missing=True,
        #     ),
        #     target_region=regions["target"],
        # ),
        # expand(
        #     expand(
        #         "results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
        #         zip,
        #         ID=samples["ID"],
        #         GENOME=samples["genome_build"],
        #         allow_missing=True,
        #     ),
        #     target_region=regions["target"],
        # ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),

        # #add length
        # expand(expand(
        #     "results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank0.csv.gz",
        #     zip,
        #     SAMPLE=samples["sample"],
        #     ID=samples["ID"],
        #     GENOME=samples["genome_build"],
        #     allow_missing=True,
        # ),target_region=regions["target"],
        # ),
        # expand(expand(
        #     "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_length.csv.gz",
        #     zip,
        #     SAMPLE=samples["sample"],
        #     ID=samples["ID"],
        #     GENOME=samples["genome_build"],
        #     allow_missing=True,
        # ),target_region=regions["target"],
        # ),

        # ## add normalize
        # expand(expand(
        #     "results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_{type}_{normalize}.tsv.gz",
        #     zip,
        #     SAMPLE=samples["sample"],
        #     ID=samples["ID"],
        #     GENOME=samples["genome_build"],
        #     allow_missing=True,
        # ),target_region=regions["target"],type=["WPS","COV","STARTS"],normalize=["normalizedFlank"] # "normalized",
        # ),

        #expand(
        #    expand(
        #        "results/plots/overlays/{ID}/{target_region}--{SAMPLE}_overlays.pdf",
        #        zip,
        #        SAMPLE=samples["sample"],
        #        ID=samples["ID"],
        #        allow_missing=True,
        #    ),
        #    target_region=regions["target"],
        #),

# genomic bam: saturated calculation, 5 fold extend to both sides, as whole transcript (error still exist if extended region exceeded boundary)
# rule add_flanks:
#     input:
#         Region_file=lambda wildcards: regions["path"][wildcards.target_region],
#     output:
#         "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
#     params:
#         #flank=1000, # orignail fixed version
#         # flank=lambda wildcards, input: get_length(input.Region_file)*3, # flank length: 10 folds of target length
#         chrsize = config['genome']['chrsize'], #config[wildcards.GENOME]['chrsize']
#         tmp1 = "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.tmp"
#     shell:
#         """
#         (fname={input.Region_file};
#         if [[ $fname == *.gz ]];
#         then zcat $fname;
#         else cat $fname;
#         fi;) > {params.tmp1}; 
#         bedtools slop -i {params.tmp1} -g {params.chrsize} -b 5 -pct | \
#         gzip -c > {output}
#         rm -rf {params.tmp1}
#         """


# tx bam: target region not exceeded (ez to exceed tx boundary), full-length tx shuf 5 regions as bg (can overlap with each other, as piRNA tx are too short)
rule get_backgroundFlank:
    input:
        region=lambda wildcards: regions["path"][wildcards.target_region],
        genome=config["genome"]["chrsize"],
        #gap=lambda wildcards: config[wildcards.GENOME]["universal_blacklist"],
    output:
        o0="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank0_regions.bed.gz",
        o1="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank1_regions.bed.gz",
        o2="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank2_regions.bed.gz",
        o3="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank3_regions.bed.gz",
        o4="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank4_regions.bed.gz",
    params:
        #length=lambda wildcards, input: get_length(input.region),
        insideExtRegion="" # if not config['tx_mapping'] else "-incl"
    conda:
        "envs/background.yml"
    shell:
        """
        bedtools shuffle {params.insideExtRegion} -seed 1 -chrom -i {input.region} -g {input.genome} | \
        gzip -c > {output.o0}
        bedtools shuffle {params.insideExtRegion} -seed 2 -chrom -i {input.region} -g {input.genome} | \
        gzip -c > {output.o1}
        bedtools shuffle {params.insideExtRegion} -seed 3 -chrom -i {input.region} -g {input.genome} | \
        gzip -c > {output.o2}
        bedtools shuffle {params.insideExtRegion} -seed 4 -chrom -i {input.region} -g {input.genome} | \
        gzip -c > {output.o3}
        bedtools shuffle {params.insideExtRegion} -seed 5 -chrom -i {input.region} -g {input.genome} | \
        gzip -c > {output.o4}
        """


# rule get_backgroundFlank:
#     input:
#         Region_file=lambda wildcards: regions["path"][wildcards.target_region],
#     output:
#         up="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank1_regions.bed.gz",
#         dn="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank2_regions.bed.gz",
#     params:
#         #flank = 1000, # orignail fixed version
#         flank = lambda wildcards, input: get_length(input.Region_file)*3, # upstream/left flank boundary: 10 folds to 5 folds of target length
#         chrsize = config['genome']['chrsize'], #config[wildcards.GENOME]['chrsize']
#     shell:
#         """        
#         (fname={input.Region_file};
#         if [[ $fname == *.gz ]];
#         then zcat $fname;
#         else cat $fname;
#         fi;) | \
#         bedtools shift -i stdin -g {params.chrsize} -s -{params.flank} | \
#         gzip -c > {output.up}

#         (fname={input.Region_file};
#         if [[ $fname == *.gz ]];
#         then zcat $fname;
#         else cat $fname;
#         fi;) | \
#         bedtools shift -i stdin -g {params.chrsize} -s +{params.flank} | \
#         gzip -c > {output.dn}
#         """


rule extract_counts:
    input:
        # target="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        target=lambda wildcards: regions["path"][wildcards.target_region],
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
        WPS_v2="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS_v2.csv.gz",
        COV="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
        STARTS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
        length="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_length.csv.gz",
    log:
        "results/intermediate/{ID}/table/{GENOME}/target/log/{target_region}--{SAMPLE}.log"
    params:
        minRL=config["minRL"], # usually target len + protection len (1bp TSS + 120bp = 121bp)
        maxRL=config["maxRL"],
        libType=lambda wildcards: samples["libType"][wildcards.SAMPLE],
        bpProtection=config["bpProtection"],
        downSampleRatio=config["downSampleRatio"],
        out_pre="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_%s.csv.gz",
    conda:
        "envs/cfDNA.yml"
    shell:
        """
        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        --protection={params.bpProtection} \
        --libType={params.libType} \
        --downsample={params.downSampleRatio} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE} > {log} 2>&1
        """


rule extract_counts_backgroundFlank:
    input:
        background0="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank0_regions.bed.gz",
        background1="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank1_regions.bed.gz",
        background2="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank2_regions.bed.gz",
        background3="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank3_regions.bed.gz",
        background4="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_backgroundFlank4_regions.bed.gz",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank0.csv.gz",
        WPS0_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank0.csv.gz",
        COV0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank0.csv.gz",
        STARTS0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank0.csv.gz",
        length0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank0.csv.gz",

        WPS1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank1.csv.gz",
        WPS1_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank1.csv.gz",
        COV1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank1.csv.gz",
        STARTS1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank1.csv.gz",
        length1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank1.csv.gz",

        WPS2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank2.csv.gz",
        WPS2_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank2.csv.gz",
        COV2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank2.csv.gz",
        STARTS2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank2.csv.gz",
        length2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank2.csv.gz",

        WPS3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank3.csv.gz",
        WPS3_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank3.csv.gz",
        COV3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank3.csv.gz",
        STARTS3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank3.csv.gz",
        length3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank3.csv.gz",

        WPS4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank4.csv.gz",
        WPS4_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank4.csv.gz",
        COV4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank4.csv.gz",
        STARTS4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank4.csv.gz",
        length4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.backgroundFlank4.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        libType=lambda wildcards: samples["libType"][wildcards.SAMPLE],
        bpProtection=config["bpProtection"],
        downSampleRatio=config["downSampleRatio"],
        out_pre0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.backgroundFlank0.csv.gz",
        out_pre1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.backgroundFlank1.csv.gz",
        out_pre2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.backgroundFlank2.csv.gz",
        out_pre3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.backgroundFlank3.csv.gz",
        out_pre4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.backgroundFlank4.csv.gz",
    log:
        log0 = "results/intermediate/{ID}/table/{GENOME}/background/log/{target_region}--{SAMPLE}.0.log",
        log1 = "results/intermediate/{ID}/table/{GENOME}/background/log/{target_region}--{SAMPLE}.1.log",
        log2 = "results/intermediate/{ID}/table/{GENOME}/background/log/{target_region}--{SAMPLE}.2.log",
        log3 = "results/intermediate/{ID}/table/{GENOME}/background/log/{target_region}--{SAMPLE}.3.log",
        log4 = "results/intermediate/{ID}/table/{GENOME}/background/log/{target_region}--{SAMPLE}.4.log",
    conda:
        "envs/cfDNA.yml"
    shell:
        """
        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
            --minInsert={params.minRL} \
            --maxInsert={params.maxRL} \
            --protection={params.bpProtection} \
            --libType={params.libType} \
            --downsample={params.downSampleRatio} \
            -i {input.background0} \
            -o {params.out_pre0} {input.BAMFILE} > {log.log0} 2>&1 

        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
            --minInsert={params.minRL} \
            --maxInsert={params.maxRL} \
            --protection={params.bpProtection} \
            --libType={params.libType} \
            --downsample={params.downSampleRatio} \
            -i {input.background1} \
            -o {params.out_pre1} {input.BAMFILE} > {log.log1} 2>&1 

        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
            --minInsert={params.minRL} \
            --maxInsert={params.maxRL} \
            --protection={params.bpProtection} \
            --libType={params.libType} \
            --downsample={params.downSampleRatio} \
            -i {input.background2} \
            -o {params.out_pre2} {input.BAMFILE} > {log.log2} 2>&1 

        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
            --minInsert={params.minRL} \
            --maxInsert={params.maxRL} \
            --protection={params.bpProtection} \
            --libType={params.libType} \
            --downsample={params.downSampleRatio} \
            -i {input.background3} \
            -o {params.out_pre3} {input.BAMFILE} > {log.log3} 2>&1 

        scripts/WPS/extractFromBAM_RegionBed_WPS_RBP_Cov.py \
            --minInsert={params.minRL} \
            --maxInsert={params.maxRL} \
            --protection={params.bpProtection} \
            --libType={params.libType} \
            --downsample={params.downSampleRatio} \
            -i {input.background4} \
            -o {params.out_pre4} {input.BAMFILE} > {log.log4} 2>&1 
        """

## add normalize by target region specifc flank regions
#results/intermediate/lulab/table/GRCh38/background/Exon1end--STAD-PKU-9-wgs_COV.background.csv.gz
rule normalize_WPS_by_flankregions:
    input:  
        target_WPS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
        target_WPS_v2="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS_v2.csv.gz",
        target_COV="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
        target_STARTS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",

        background_WPS0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank0.csv.gz",
        background_WPS0_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank0.csv.gz",
        background_COV0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank0.csv.gz",
        background_STARTS0="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank0.csv.gz",

        background_WPS1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank1.csv.gz",
        background_WPS1_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank1.csv.gz",
        background_COV1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank1.csv.gz",
        background_STARTS1="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank1.csv.gz",

        background_WPS2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank2.csv.gz",
        background_WPS2_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank2.csv.gz",
        background_COV2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank2.csv.gz",
        background_STARTS2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank2.csv.gz",

        background_WPS3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank3.csv.gz",
        background_WPS3_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank3.csv.gz",
        background_COV3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank3.csv.gz",
        background_STARTS3="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank3.csv.gz",

        background_WPS4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.backgroundFlank4.csv.gz",
        background_WPS4_v2="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS_v2.backgroundFlank4.csv.gz",
        background_COV4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.backgroundFlank4.csv.gz",
        background_STARTS4="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.backgroundFlank4.csv.gz",
    output:
        output_WPS="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_WPS_normalizedFlank.tsv.gz",
        output_WPS_v2="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_WPS_v2_normalizedFlank.tsv.gz",
        output_WPS_v3="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_WPS_v3_normalizedFlank.tsv.gz",
        output_COV="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_COV_normalizedFlank.tsv.gz",
        output_STARTS="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_STARTS_normalizedFlank.tsv.gz",
    log:
        # log="results/intermediate/{ID}/table/{GENOME}/normalize/log/{target_region}--{SAMPLE}.log"
    conda:
        "envs/overlays.yml"
    script:
        """scripts/expression_analysis/normalize_byflank_tx.py"""
