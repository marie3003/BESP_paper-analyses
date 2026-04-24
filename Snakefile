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
def get_summary_trees(wildcards):
    import pandas as pd
    csv = checkpoints.check_mcmc.get().output.csv
    df = pd.read_csv(csv)
    return [p.replace(".trees", "_summary.tree") for p in df["trees_path"]]


rule all:
    input:
        get_summary_trees


rule all_alignments:
    input:
        expand(
            f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}_{{i}}_snps.fasta",
            sampling=SAMPLING_TYPES, popmodel=POP_MODELS, i=range(NREPLICATES)
        ),
        f"{SIMDATA}/snp_summary.csv"


rule all_beast_test:
    input:
        expand(
            f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/seed{{seed}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T0.log",
            model=INFERENCE, sampling=SAMPLING_TYPES, popmodel=POP_MODELS, mutsig=MUTSIGS, seed=SEEDS
        )


rule all_xmls:
    input:
        expand(
            f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.xml",
            model=INFERENCE, sampling=SAMPLING_TYPES, popmodel=POP_MODELS, mutsig=MUTSIGS, i=range(NREPLICATES)
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
    log:
        "logs/simulate_trees/simulate_trees.log"
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        "r/4.4.0",
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 300,     # 5h
        cpus_per_task = 16,
    shell:
        "Rscript {input.script} > {log} 2>&1"


# =============================================================================
# Step 2a: Simulate alignment with seq-gen and extract SNPs with snp-sites
#         Full fasta is deleted immediately after SNP extraction to save space
# =============================================================================
rule simulate_alignment:
    input:
        trees = f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}.trees",
    output:
        snps = f"{SIMDATA}/{{sampling}}/{{popmodel}}/{{popmodel}}_{{i}}_snps.fasta",
    log:
        "logs/simulate_alignment/{sampling}_{popmodel}_{i}.log"
    resources:
        runtime = lambda wildcards: 60 if wildcards.sampling == "independenthomochronous" else 180,
        cpus_per_task = 1,
        mem_mb_per_cpu = 8000,
    shell:
        """
        TREE=$(sed -n "$(( {wildcards.i} + 1 ))p" {input.trees})
        TMPNWK=$(mktemp --suffix=.nwk)
        TMPFASTA=$(mktemp --suffix=.fasta)
        echo -e "$TREE\n" > $TMPNWK
        {SEQGEN} -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l10000000 -s4.6e-8 -n1 \
            < $TMPNWK > $TMPFASTA 2>> {log}
        rm $TMPNWK
        conda run -n snp_sites snp-sites -o {output.snps} $TMPFASTA >> {log} 2>&1
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
    log:
        "logs/snp_summary/snp_summary.log"
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
# Step 3: Generate BEAST XML files
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
        expand(
            f"{BEAST_DIR}/{{{{model}}}}/{{{{sampling}}}}/{{{{popmodel}}}}/{{{{mutsig}}}}/{{{{model}}}}_{{{{sampling}}}}_{{{{popmodel}}}}_{{{{mutsig}}}}.T{{i}}.xml",
            i=range(NREPLICATES)
        )
    log:
        "logs/make_beast_xml/{model}_{sampling}_{popmodel}_{mutsig}.log"
    resources:
        mem_mb_per_cpu = 2000,
        runtime = 5, # 5 min per scenario (100 xml files per scenario)
        cpus_per_task = 1,
    shell:
        "conda run -n beast_tools python scripts/MakeBEASTXML.py -c {input.config} > {log} 2>&1"


# =============================================================================
# Step 4: Run BEAST (one job per XML file per seed)
# =============================================================================
rule run_beast:
    input:
        xml = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.xml",
    output:
        log   = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/seed{{seed}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.log",
        trees = f"{BEAST_DIR}/{{model}}/{{sampling}}/{{popmodel}}/{{mutsig}}/seed{{seed}}/{{model}}_{{sampling}}_{{popmodel}}_{{mutsig}}.T{{i}}.trees",
    log:
        "logs/run_beast/{model}_{sampling}_{popmodel}_{mutsig}_T{i}_seed{seed}.log"
    benchmark:
        "benchmarks/run_beast/{model}_{sampling}_{popmodel}_{mutsig}_T{i}_seed{seed}.tsv"
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle/3.1.2",
    resources:
        mem_mb_per_cpu = 1000,
        runtime = lambda wildcards, attempt: (480 if wildcards.model == "constcoal" else 960) * (2 if wildcards.mutsig == "highmutsig" else 1),
        cpus_per_task = 2,
        slurm_jobname = lambda wildcards: f"{wildcards.model}_{wildcards.sampling}_{wildcards.popmodel}_{wildcards.mutsig}_T{wildcards.i}_s{wildcards.seed}",
    shell:
        """
        mkdir -p $(dirname {output.log})
        cd $(dirname {output.log}) && beast -overwrite -seed {wildcards.seed} $OLDPWD/{input.xml} > $OLDPWD/{log} 2>&1
        """


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
    log:
        "logs/combine_runs/{model}_{sampling}_{popmodel}_{mutsig}_T{i}.log"
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
        "libbeagle/3.1.2",
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 5,
        cpus_per_task = 1,
    shell:
        """
        mkdir -p $(dirname {output.log})
        logcombiner -burnin 1000 {input.logs} {output.log} > {log} 2>&1
        logcombiner -trees -burnin 1000 {input.trees} {output.trees} >> {log} 2>&1
        """


# =============================================================================
# Step 6: Check MCMC convergence and create summary trees.
# Checkpoint is needed here because we first need to know which runs were
# successful before we can create the summary trees in the next step.
# In all other cases we knew which files to expect from the start.
# =============================================================================
checkpoint check_mcmc:
    input:
        script = "scripts/evaluate_mcmc.R",
    output:
        csv = "scripts/successful_mcmc_runs.csv",
    log:
        "logs/check_mcmc/check_mcmc.log"
    envmodules:
        "stack/2024-06",
        "r/4.4.0",
    resources:
        mem_mb_per_cpu = 2000,
        runtime = 60,
        cpus_per_task = 16,
    shell:
        "Rscript {input.script} > {log} 2>&1"


# =============================================================================
# Step 7: Run TreeAnnotator on each successful run
# =============================================================================
rule make_summary_tree:
    input:
        trees = "{path}.trees",
        csv   = "scripts/successful_mcmc_runs.csv",
    output:
        "{path}_summary.tree",
    log:
        "logs/make_summary_tree/{path}_summary.log"
    envmodules:
        "stack/2024-06",
        "openjdk/21.0.3_9",
        "gcc/12.2.0",
        "beast1/1.10.4",
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 30,
        cpus_per_task = 1,
    shell:
        "treeannotator -burnin 0 -heights mean {input.trees} {output} > {log} 2>&1"
