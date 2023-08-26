<!-- vscode-markdown-toc -->
* 1. [ Usage](#Usage)
	* 1.1. [Step 1: Obtain a copy of this workflow](#Step1:Obtainacopyofthisworkflow)
	* 1.2. [Step 2: Configure workflow](#Step2:Configureworkflow)
	* 1.3. [Step 3: Install Snakemake](#Step3:InstallSnakemake)
	* 1.4. [Step 4: Execute workflow](#Step4:Executeworkflow)
* 2. [check output dir](#checkoutputdir)
* 3. [Acknowledgements](#Acknowledgements)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc --># Snakemake workflow: Analysis of cfRNA signals captured by fragmentation patterns of cell-free RNA



##  1. <a name='Usage'></a> Usage

###  1.1. <a name='Step1:Obtainacopyofthisworkflow'></a>Step 1: Obtain a copy of this workflow

```bash
#Clone the repository to the place to perform the data analysis
git clone https://github.com/hunnngry/tWPS.git
```

###  1.2. <a name='Step2:Configureworkflow'></a>Step 2: Configure workflow
```bash
dst="GSE71008_NCpool"

#specify required parameters and file paths
vi config/$dst/config.yml

#specify your sample setup
vi config/$dst/samples.tsv

#specify target regions
vi config/$dst/regions.tsv
```


###  1.3. <a name='Step3:InstallSnakemake'></a>Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -c bioconda -c conda-forge -n tWPS snakemake
```


###  1.4. <a name='Step4:Executeworkflow'></a>Step 4: Execute workflow

The workflows are executed from the repository root folder. The different analyses have to be executed separately. To specify the respective workflow use the `-s` switch followed by the path of the `Snakefile` (e.g.: `./snakemake/test.smk`)


```bash
#Activate the conda environment
conda activate tWPS
mkdir -p log/${dst}

#Test your configuration by performing a dry-run via
snakemake \
    --snakefile snakemake/snakefile_WPS_RBP_tx.smk \
    --configfile config/${dst}/config.yml \
    -npr 

#op1: Run in local machine (PC/MAC)
snakemake --rerun-incomplete --keep-going --printshellcmds --reason --use-conda --nolock --latency-wait 20 --restart-times 1 --jobs 6  \
    --snakefile snakemake/snakefile_WPS_RBP_tx.smk \
    --configfile config/${dst}/config.yml \
    > log/${dst}/run-${dst}_local.log 2>&1

#op2: Run in cluster (bsub)
snakemake --rerun-incomplete --keep-going --printshellcmds --reason --use-conda --nolock --latency-wait 80 --restart-times 1 --jobs 100  \
    --snakefile snakemake/snakefile_WPS_RBP_tx.smk \
    --configfile config/${dst}/config.yml \
    --cluster-config snakemake/cluster/cluster-lsf.json \
    --cluster "bsub -n {cluster.threads} -J {cluster.jobname} -q {cluster.queue} -o {cluster.output} -e {cluster.error} " \
    > log/${dst}/run-${dst}_cluster.log 2>&1
```

##  2. <a name='checkoutputdir'></a>check output dir
```bash
ls -hl results/intermediate/$dst
```

##  3. <a name='Acknowledgements'></a>Acknowledgements
tWPS for cfRNA was adapted from WPS for cfDNA [repository](https://github.com/shendurelab/cfDNA)