# =============================================================================
# Snakefile for BESP simulation study
# Run from project root: BESP_paper-analyses/
#
# Usage:
#   Dry run:   snakemake -n
#   Cluster:   snakemake --profile profiles/euler --rerun-incomplete
# =============================================================================

SAMPLING_TYPES = ["independenthomochronous", "linearconstant"]
POP_MODELS     = ["expgrowthfast", "expgrowthslow", "uniform", "bottleneck"]
NREPLICATES    = 100
INFERENCE      = ["constcoal", "skyline"]
MUTSIGS        = ["lowmutsig", "medmutsig", "highmutsig"]
SEEDS          = [44, 45, 46]   # one BEAST run per seed, combined later

SIMDATA    = "results/run1/simulated_data"
BEAST_DIR  = "results/run1/beast_inference"
CONFIG_DIR = "results/run1/config"
SEQGEN     = workflow.basedir + "/../Seq-Gen-1.3.5/source/seq-gen"

def config_file(wildcards):
    return f"{CONFIG_DIR}/{wildcards.model}_{wildcards.sampling}_{wildcards.popmodel}_{wildcards.mutsig}.cfg"


# =============================================================================
# Rule all — final targets
# =============================================================================
rule all:
    input:
        # Summary trees for all successful MCMC runs (end of pipeline)
        "scripts/successful_mcmc_runs.csv"


rule all_alignments:
    input:
        expand(
            f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}_{{i}}_snps.fasta",
            sampling=SAMPLING_TYPES, popmodel=POP_MODELS, i=range(NREPLICATES)
        ),
        f"{SIMDATA}/snp_summary.csv"


rule all_xmls:
    input:
        expand(
            f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}",
            model=INFERENCE, sampling=SAMPLING_TYPES, popmodel=POP_MODELS, mutsig=MUTSIGS
        )


# =============================================================================
# Step 1: Simulate trees
# =============================================================================
rule simulate_trees:
    input:
        script = "scripts/simulate_trees.R",
        utils  = "scripts/SimUtils.R",
    output:
        expand(
            f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}.trees",
            sampling=SAMPLING_TYPES, popmodel=POP_MODELS
        )
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        "r/4.4.0",
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 300,     # 5h
        cpus_per_task = 16,
    shell:
        "Rscript {input.script}"


# =============================================================================
# Step 2a: Simulate alignment with seq-gen and extract SNPs with snp-sites
#         Full fasta is deleted immediately after SNP extraction to save space
# =============================================================================
rule simulate_alignment:
    input:
        trees = f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}.trees",
    output:
        snps = f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}_{{i}}_snps.fasta",
    resources:
        runtime = 180,    # 3h
        cpus_per_task = 1,
        mem_mb_per_cpu = 8000,
    shell:
        """
        TREE=$(sed -n "$(( {wildcards.i} + 1 ))p" {input.trees})
        TMPNWK=$(mktemp --suffix=.nwk)
        TMPFASTA=$(mktemp --suffix=.fasta)
        echo -e "$TREE\n" > $TMPNWK
        {SEQGEN} -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l10000000 -s4.6e-8 -n1 \
            < $TMPNWK > $TMPFASTA
        rm $TMPNWK
        conda run -n snp_sites snp-sites -o {output.snps} $TMPFASTA
        rm $TMPFASTA
        """

# =============================================================================
# Step 2b: Summarize number of SNPs in each alignment and save to CSV
# =============================================================================
rule snp_summary:
    input:
        expand(
            f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}_{{i}}_snps.fasta",
            sampling=SAMPLING_TYPES, popmodel=POP_MODELS, i=range(NREPLICATES)
        )
    output:
        f"{SIMDATA}/snp_summary.csv"
    resources:
        mem_mb_per_cpu = 2000,
        runtime = 10,
        cpus_per_task = 1,
    run:
        import os
        rows = []
        for f in input:
            parts = f.split("/")
            sampling = parts[-3]
            popmodel = parts[-2]
            i = os.path.basename(f).replace("_snps.fasta", "").split("_")[-1]
            with open(f) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        n_snps = len(line.strip())
                        break
            rows.append(f"{sampling},{popmodel},{i},{n_snps}")
        with open(output[0], "w") as out:
            out.write("sampling,popmodel,replicate,n_snps\n")
            out.write("\n".join(rows) + "\n")


# =============================================================================
# Step 3: Generate BEAST XML files for one scenario
# =============================================================================
rule make_beast_xml:
    input:
        snps   = expand(
                     f"{SIMDATA}/{{{{sampling}}}}/{{{{popmodel}}}}/{{{{popmodel}}}}_{{i}}_snps.fasta",
                     i=range(NREPLICATES)
                 ),
        trees  = f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}.trees",
        config = config_file,
    output:
        directory(f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}")
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 30,
        cpus_per_task = 1,
    shell:
        "conda run -n beast_tools python scripts/MakeBEASTXML.py -c {input.config}"


# =============================================================================
# Step 4: Run BEAST (one job per XML file per seed)
# =============================================================================
rule run_beast:
    input:
        xml     = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.xml",
        xml_dir = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}"
    output:
        log   = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/seed{{seed}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.log",
        trees = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/seed{{seed}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.trees",
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle",
    resources:
        mem_mb_per_cpu = 8000,
        runtime = 960,   # 16h
        cpus_per_task = 1,
    shell:
        "beast -overwrite -seed {wildcards.seed} -working {input.xml}"


# =============================================================================
# Step 5: Combine BEAST runs with logcombiner
# =============================================================================
rule combine_runs:
    input:
        logs  = lambda wc: expand(
                    f"{BEAST_DIR}/{wc.model}/{wc.sampling}/{wc.popmodel}/{wc.mutsig}/seed{{seed}}/{wc.model}_{wc.sampling}_{wc.popmodel}_{wc.mutsig}.T{wc.i}.log",
                    seed=SEEDS
                ),
        trees = lambda wc: expand(
                    f"{BEAST_DIR}/{wc.model}/{wc.sampling}/{wc.popmodel}/{wc.mutsig}/seed{{seed}}/{wc.model}_{wc.sampling}_{wc.popmodel}_{wc.mutsig}.T{wc.i}.trees",
                    seed=SEEDS
                ),
    output:
        log   = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/combined/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.log",
        trees = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/combined/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.trees",
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle",
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 10,
        cpus_per_task = 1,
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
        mem_mb_per_cpu = 2000,
        runtime = 60,
        cpus_per_task = 16,
    shell:
        "Rscript {input.script}"
