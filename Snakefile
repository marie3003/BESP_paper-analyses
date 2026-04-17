# =============================================================================
# Snakefile for BESP simulation study
# Run from project root: BESP_paper-analyses/
#
# Usage:
#   Dry run:   snakemake -n
#   Cluster:   snakemake --profile profiles/euler --rerun-incomplete
# =============================================================================

SAMPLING_TYPES = ["independent_homochronous", "linear_constant"]
POP_MODELS     = ["expgrowth_fast", "expgrowth_slow", "uniform", "bottleneck"]
NREPLICATES    = 100
INFERENCE      = ["constcoal", "skyline"]
MUTSIGS        = ["lowmutsig", "medmutsig", "highmutsig"]
SEEDS          = [44, 45, 46]   # one BEAST run per seed, combined later

SIMDATA    = "results/run1/simulated_data"
BEAST_DIR  = "results/run1/beast_inference"
CONFIG_DIR = "results/run1/config"
SEQGEN     = "Seq-Gen-1.3.5/source/seq-gen"

def short(name):
    return name.replace("_", "")

def config_file(wildcards):
    return (f"{CONFIG_DIR}/{wildcards.model}_{short(wildcards.sampling)}"
            f"_{short(wildcards.pop)}_{wildcards.mutsig}.cfg")

def beast_log(wildcards):
    return (f"{BEAST_DIR}/{wildcards.model}/{wildcards.sampling}/{wildcards.pop}"
            f"/{wildcards.mutsig}/{wildcards.model}_{short(wildcards.sampling)}"
            f"_{short(wildcards.pop)}_{wildcards.mutsig}.T{wildcards.i}.log")


# =============================================================================
# Rule all — final targets
# =============================================================================
rule all:
    input:
        # Summary trees for all successful MCMC runs (end of pipeline)
        "scripts/successful_mcmc_runs.csv"


# =============================================================================
# Step 1: Simulate trees
# =============================================================================
rule simulate_trees:
    input:
        script = "scripts/simulate_trees.R",
        utils  = "scripts/SimUtils.R",
    output:
        expand(
            f"{SIMDATA}/{{sampling}}/{{pop}}/{{pop}}.trees",
            sampling=SAMPLING_TYPES, pop=POP_MODELS
        )
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        "r/4.4.0",
    resources:
        mem_mb  = 64000,   # 16 CPUs x 4G
        runtime = 300,     # 5h
        cpus    = 16,
    shell:
        "Rscript {input.script}"


# =============================================================================
# Step 2: Simulate alignment with seq-gen and extract SNPs with snp-sites
#         Full fasta is deleted immediately after SNP extraction to save space
# =============================================================================
rule simulate_alignment:
    input:
        trees = f"{SIMDATA}/{{sampling}}/{{pop}}/{{pop}}.trees",
    output:
        snps = f"{SIMDATA}/{{sampling}}/{{pop}}/{{pop}}_{{i}}_snps.fasta",
    resources:
        mem_mb  = 8000,
        runtime = 180,    # 3h
        cpus    = 1,
    shell:
        """
        TREE=$(sed -n "$(( {wildcards.i} + 1 ))p" {input.trees})
        TMPNWK=tmp_tree_{wildcards.sampling}_{wildcards.pop}_{wildcards.i}.nwk
        TMPFASTA=tmp_fasta_{wildcards.sampling}_{wildcards.pop}_{wildcards.i}.fasta
        echo -e "$TREE\n" > $TMPNWK
        {SEQGEN} -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l10000000 -s4.6e-8 -n1 \
            < $TMPNWK > $TMPFASTA
        rm $TMPNWK
        conda run -n snp_sites snp-sites -o {output.snps} $TMPFASTA
        rm $TMPFASTA
        """


# =============================================================================
# Step 3: Generate BEAST XML files for one scenario
# =============================================================================
rule make_beast_xml:
    input:
        snps   = expand(
                     f"{SIMDATA}/{{{{sampling}}}}/{{{{pop}}}}/{{{{pop}}}}_{{i}}_snps.fasta",
                     i=range(NREPLICATES)
                 ),
        trees  = f"{SIMDATA}/{{sampling}}/{{pop}}/{{pop}}.trees",
        config = config_file,
    output:
        directory(f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}")
    resources:
        mem_mb  = 4000,
        runtime = 30,
        cpus    = 1,
    shell:
        "conda run -n beast_tools python scripts/MakeBEASTXML.py -c {input.config}"


# =============================================================================
# Step 4: Run BEAST (one job per XML file per seed)
# =============================================================================
rule run_beast:
    input:
        xml = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}/{{model}}_{short('{{sampling}}')}_{{pop_short}}_{{mutsig}}.T{{i}}.xml",
        xml_dir = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}"
    output:
        log   = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}/seed{{seed}}/{{model}}_{short('{{sampling}}')}_{{pop_short}}_{{mutsig}}.T{{i}}.log",
        trees = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}/seed{{seed}}/{{model}}_{short('{{sampling}}')}_{{pop_short}}_{{mutsig}}.T{{i}}.trees",
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle",
    resources:
        mem_mb  = 8000,
        runtime = 960,   # 16h
        cpus    = 1,
    shell:
        "beast -overwrite -seed {wildcards.seed} -working {input.xml}"


# =============================================================================
# Step 5: Combine BEAST runs with logcombiner
# =============================================================================
rule combine_runs:
    input:
        logs  = expand(
                    f"{BEAST_DIR}/{{{{model}}}}/{{{{sampling}}}}/{{{{pop}}}}/{{{{mutsig}}}}/seed{{seed}}/{{{{model}}}}_{{short('{{{{sampling}}}}')}}_{{}}_{{{{mutsig}}}}.T{{{{i}}}}.log",
                    seed=SEEDS
                ),
        trees = expand(
                    f"{BEAST_DIR}/{{{{model}}}}/{{{{sampling}}}}/{{{{pop}}}}/{{{{mutsig}}}}/seed{{seed}}/{{{{model}}}}_{{short('{{{{sampling}}}}')}}_{{}}_{{{{mutsig}}}}.T{{{{i}}}}.trees",
                    seed=SEEDS
                ),
    output:
        log   = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}/combined/{{name}}.T{{i}}.log",
        trees = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{pop}}/{{mutsig}}/combined/{{name}}.T{{i}}.trees",
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle",
    resources:
        mem_mb  = 4000,
        runtime = 10,
        cpus    = 1,
    shell:
        """
        mkdir -p $(dirname {output.log})
        logcombiner -burnin 1000 {input.logs} {output.log}
        logcombiner -trees -burnin 1000 {input.trees} {output.trees}
        """


# =============================================================================
# Step 6: Check MCMC convergence and create summary trees
# =============================================================================
rule check_mcmc:
    input:
        script = "scripts/evaluate_mcmc.R",
    output:
        csv = "scripts/successful_mcmc_runs.csv",
    envmodules:
        "stack/2024-06",
        "r/4.4.0",
    resources:
        mem_mb  = 32000,   # 16 CPUs x 2G
        runtime = 60,
        cpus    = 16,
    shell:
        "Rscript {input.script}"
