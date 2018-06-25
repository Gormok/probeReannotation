### As example, the pipeline is currently ran on a select subset of probes
### To rebuild all, use -F snakemake flag

### REQUIRED ###
1. Anaconda : https://anaconda.org
2. Human Genome (GTF, FA)
3. Index files for tools (BOWTIE, STAR)

### FOLDER STRUCTURE ###
Input folders:
    Raw: Probe descriptions (manufacturer)
    Orig: Original symbol lists (manufacturer)
    selectProbes: Probe subset - for testing

Computation folders:
    ProbeSequences: Fastq files derived from Raw
    Alignment: Probe alignment results (sam)
    SymbolLists: Reannotated probe symbol lists
    Comparisons: Comparisons of reannotations
        mismatches: mismatched annotations (w.r.t. Orig)
        refmatch: matched annotations (w.r.t. Orig)
        total: total annotated probe counts
        genes: total unique ENSEMBL ids covered
        total-refmatch: 2D plot comparing the two values


### USAGE ###
1. Fill in paths on your system for genome and index paths in file Setup.yaml
    GTF, FA, BOWTIE, STAR

2. Create the base environment 'probeAlignment'
    'conda env create --file envs/env1.yaml'
3. Run the pipeline
    single core: 'snakemake --snakefile Comparisons.rules --use-conda'
    multiple cores (8): 'snakemake --snakefile Comparisons.rules --use-conda --cores 8'
    [ '--use-conda' required for tophat dependency (envs/env2.yaml) ]




The DAG ('dag.svg') can be created by
    'snakemake -s Comparisons.rules --dag | dot -Tsvg > dag.svg'
