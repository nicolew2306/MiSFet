# vim: syntax=python expandtab ts=4
# Snakemake for MIDAS2-Signature-FEAST

# Config file for the values
configfile: "sm_config/config.yaml"
INPUTDIR = Path(config["inputdir"])


# Adjust the glob pattern to fit your input file names
SAMPLES = glob_wildcards("input/{sample}_1.fq.gz").sample
#SAMPLES = set(glob_wildcards(INPUTDIR/config["input_fn_pattern"]).sample)
INPUT_read1 = INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="1"),
INPUT_read2 = INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="2")
print("Found the following samples:")
print(SAMPLES)
species_list_file=config["species_list_file"]
SPECIES_IDS=[sid.strip() for sid in open(species_list_file).readlines()]

    
rule init_db:
    """
    initialize midas2 gtdb database
    """
    output:
        db="my_midasdb_gtdb/metadata.tsv",
    log:
        stderr="logs/init_db.stderr",
    threads:
        2
    conda:
        "envs/midas2.yml"
    shell:
        """
        midas2 database \
            --init \
            --midasdb_name gtdb \
            --midasdb_dir my_midasdb_gtdb \
            2> {log.stderr} 
        """


rule run_species:
        """
        Using the basic markers to run species from the reads.
        """
        input:
            read1="input/{sample}_1.fq.gz",
            read2="input/{sample}_2.fq.gz",
            #read1=INPUT_read1,
            #read2=INPUT_read2,
            db="my_midasdb_gtdb/metadata.tsv",
        output:
            sp="output/{sample}/species/markers_profile.tsv",
        log:
            stderr="logs/run_species/{sample}.stderr"
        threads:
            20
        conda:
            "envs/midas2.yml"
        shell:
           """
            midas2 run_species \
                --sample_name {wildcards.sample} \
                -1 {input.read1} \
                -2 {input.read2} \
                --midasdb_name gtdb \
                --midasdb_dir my_midasdb_gtdb \
                --num_cores {threads} \
                output \
                2> {log.stderr}
            """

rule sample_list:
    """
    get list of sample names
    """
    output:
        "sample_list.tsv",
    shell:
        """
        echo "sample_name\tmidas_outdir" > {output}
        for sample in {SAMPLES}; do
            echo $sample"\t"output >> {output}
        done
        """ 

rule merge_species:
    """
    here we merge the species created above.
    """
    input:
        sample_list="sample_list.tsv",
        sp=expand("output/{sample}/species/markers_profile.tsv",sample=SAMPLES),
    output:
        merge_output="output/merge/species/species_prevalence.tsv",
    log:
        stderr="logs/merge.stderr"
    params:
        min_cov=config["MIN_COV"]
    threads:
        2
    conda:
        "envs/midas2.yml"
    shell:
        """
        set +u +o pipefail
        midas2 merge_species \
            --samples_list {input.sample_list} \
            --min_cov {params.min_cov} \
            output/merge \
            2> {log.stderr}
        """   

rule build_snp_db:
    """
    build specialized database
    """
    input:
        species_list=species_list_file,
        db="my_midasdb_gtdb/metadata.tsv"
    output:
        expand("my_midasdb_gtdb/pangenomes/{species_id}/cluster_info.txt",species_id=SPECIES_IDS),
    log:
        stderr="logs/build_snp_db.stderr"
    threads:
        8
    conda:
        "envs/midas2.yml"
    shell:
        """
        midas2 database \
            --download \
            --midasdb_name gtdb \
            --midasdb_dir my_midasdb_gtdb \
            --species_list {input.species_list} \
            2> {log.stderr}
        find my_midasdb_gtdb/pangenomes -exec touch {{}} \;
        """

rule run_snps:
        """
        Using the specialized markers to run snps for species from the reads.
        """
        input:
            #read1=INPUT_read1,
            #read2=INPUT_read2,
            read1="input/{sample}_1.fq.gz",
            read2="input/{sample}_2.fq.gz",
            db=expand("my_midasdb_gtdb/pangenomes/{species_id}/cluster_info.txt",species_id=SPECIES_IDS),
        output:
            sp="output/{sample}/snps/snps_summary.tsv",
        log:
            stderr="logs/run_snps/{sample}.stderr"
        threads:
            20
        conda:
            "envs/midas2.yml"
        params:
            select_by=config["select_by"]
        shell:
           """
            midas2 run_snps \
                --sample_name {wildcards.sample} \
                -1 {input.read1} \
                -2 {input.read2} \
                --midasdb_name gtdb \
                --midasdb_dir my_midasdb_gtdb \
                --select_by {params.select_by} \
                --select_threshold=2,0.5 \
                --num_cores {threads} \
                output \
                2> {log.stderr}
            """

rule merge_snps:
    """
    merging the snps.
    """
    input:
        sample_list="sample_list.tsv",
        sp=expand("output/{sample}/snps/snps_summary.tsv",sample=SAMPLES),
    output:
        merge_output="output/merge/snps/snps_summary.tsv",
    log:
        stderr="logs/snp_merge.stderr"
    threads:
        20
    conda:
        "envs/midas2.yml"
    shell:
        """
        midas2 merge_snps \
            --samples_list {input.sample_list} \
            --midasdb_name gtdb \
            --midasdb_dir my_midasdb_gtdb \
            --num_cores {threads} \
            output/merge \
            2> {log.stderr}
        """   

checkpoint merged_snps:
    """
    see which species merged.
    """
    input:
        snp_dir_init="output/merge/snps/snps_summary.tsv", #### get folder name, use output from merge_snp rule
        species_list=species_list_file,
    output:
        snp_merge_file="snps_merged.txt"
    threads:
        2
    shell:
        """
        directory_path=$(dirname {input.snp_dir_init})
        ls $directory_path | grep -xF -f {input.species_list} > {output.snp_merge_file} 
        """

rule lz4_to_bz2:
    """
    changing the snp outputs so Signature can use them.
    """
    input: 
        merge="snps_merged.txt",
        directory="output/merge/snps/snps_summary.tsv"
    output:
        snp_freq_output="output/merge/snps/{snp_id}/snps_ref_freq.txt.bz2",
        snp_depth_output="output/merge/snps/{snp_id}/snps_depth.txt.bz2",
    threads:
        2
    shell:
        """
        directory_path=$(dirname {input.directory})
        lz4 -dc $directory_path/{wildcards.snp_id}/{wildcards.snp_id}.snps_freqs.tsv.lz4 > $directory_path/{wildcards.snp_id}/snps_ref_freq.txt 
        bzip2 $directory_path/{wildcards.snp_id}/snps_ref_freq.txt
        lz4 -dc $directory_path/{wildcards.snp_id}/{wildcards.snp_id}.snps_depth.tsv.lz4 > $directory_path/{wildcards.snp_id}/snps_depth.txt
        bzip2 $directory_path/{wildcards.snp_id}/snps_depth.txt
        """

rule change_signature_config:
    """
    changing the config file for signature. But I need an output for this but output already exists.
    """
    output:
        "configs/config.yaml"
    threads:
        2
    shell:
        """
        echo "input_dir: 'metadata'" >  configs/config.yaml
        echo "snp_dir: 'output/merge/snps'" >> configs/config.yaml
        echo "output_dir: 'output/Signature_SNVs'" >>  configs/config.yaml
        """

rule running_signature:
    """
    Running Signature to get the SNP count.
    """
    input:
        config="configs/config.yaml",
        snp_freq_output="output/merge/snps/{snp_id}/snps_ref_freq.txt.bz2",
        snp_depth_output="output/merge/snps/{snp_id}/snps_depth.txt.bz2",
    output:
       directory("output/Signature_SNVs/{snp_id}")
    conda:
        "envs/signature.yaml"
    params:
        min_reads=config["MIN_READS"],
        start_index=config["START_INDEX"],
        end_index=config["END_INDEX"],
        config_file_path="configs/config.yaml"
    log:
        stderr="logs/signature/siganture_{snp_id}.stderr"
    threads:
        20
    shell:
        """
        python Signature-SNVs/src/signature_snvs/signature_snvs_cli.py \
            --species {wildcards.snp_id} \
            --min_reads {params.min_reads} \
            --start_index {params.start_index} \
            --end_index {params.end_index} \
            --config_file_path {params.config_file_path} \
            2> {log.stderr}
        """

rule run_feast:
    """
    Finally running FEAST using the R script.
    """
    input:
        metadata="metadata/sink_source.csv",
        sigs="output/Signature_SNVs/{snp_id}/"
    output:
        indicator="metadata/.{snp_id}.success"
    shadow:
        "shallow"
    conda:
        "envs/sig_to_FEAST.yaml"
    threads:
        2
    shell:
        """
        find {input.sigs} -exec bzip2 -dk {{}} \;
        Rscript -e 'devtools::install_github("cozygene/FEAST")'
        Rscript signature_to_feast.R \
            {input.sigs} \
            {input.metadata}
        touch {output.indicator}
        cp *txt {input.sigs} 
        """

def merged_snps_ids(wildcards):   
    with checkpoints.merged_snps.get().output[0].open() as f: 
        snps_ids = [snps_id.strip() for snps_id in f.readlines()]
    return [f"metadata/.{snp_id}.success" for snp_id in snps_ids]

rule initial_rule:
    """
    Initial rule to indicate that there is a checkpoint that needs to be satisfied.
    """
    default_target: True
    input:
        merged_snps_ids
