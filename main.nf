#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lp-wgs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/chelauk/lp-wgs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_lp_wgs_pipeline'
include { LP_WGS                           } from './workflows/lp_wgs'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_lp_wgs_pipeline'

//
// WORKFLOW: Run main lp-wgs analysis pipeline
//
workflow {
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GENOME PARAMETER VALUES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    //def ref_dict                  = params.dict ?: getGenomeAttribute('dict')
    def ref_fasta                 = params.fasta ?: getGenomeAttribute('fasta')
    def ref_fasta_fai             = params.fasta_fai ?: getGenomeAttribute('fasta_fai')
    def ref_bwa                   = params.bwa ?: getGenomeAttribute('bwa')
    def ref_centromere            = params.centromere ?: getGenomeAttribute('centromere')
    //def ref_map_wig               = params.map_wig ?: getGenomeAttribute('map_wig')
    def ref_map_wig_file          = params.map_wig_file ?: getGenomeAttribute('map_wig_file')
    def ref_gc_wig                = params.gc_wig ?: getGenomeAttribute('gc_wig')
    def ref_ichor_genome_build    = params.ichor_genome_build ?: getGenomeAttribute('ichor_genome_build')
    def ref_ichor_genome_style    = params.ichor_genome_style ?: getGenomeAttribute('ichor_genome_style')
    def ref_chr_bed               = params.mosdepth_bed ?: getGenomeAttribute('mosdepth_bed') ?: params.chr_bed ?: getGenomeAttribute('chr_bed')
    def ref_medicc_arms           = params.medicc_arms ?: getGenomeAttribute('medicc_arms')
    def ref_medicc_genes          = params.medicc_genes ?: getGenomeAttribute('medicc_genes')
    def ref_chr_arm_boundaries    = params.chr_arm_boundaries ?: getGenomeAttribute('chr_arm_boundaries')
    def ref_qdnaseq_genome        = params.qdnaseq_genome ?: getGenomeAttribute('qdnaseq_genome')
    def ref_qdnaseq_package       = params.qdnaseq_package ?: getGenomeAttribute('qdnaseq_package')
    
    ch_bin_dir = channel.fromPath("${projectDir}/bin/", type: 'dir').collect()  
    //def ref_hmmcopy_chromosomes   = params.hmmcopy_chromosomes ?: getGenomeAttribute('hmmcopy_chromosomes')

    // Initialise genome resources close to the workflow entrypoint.
    ch_fasta = ref_fasta ? channel.fromPath(ref_fasta).map { path -> [[id: path.baseName], path] }.collect() : channel.empty()
    ch_fasta_fai = ref_fasta_fai ? channel.fromPath(ref_fasta_fai).map { path -> [[id: path.name.replaceFirst(/\.(fa|fasta)\.fai$/, '')], path] }.collect() : channel.empty()
    //ch_dict = ref_dict ? channel.fromPath(ref_dict).collect() : channel.empty()
    ch_chr_arm_boundaries = ref_chr_arm_boundaries ? channel.fromPath(ref_chr_arm_boundaries).collect() : channel.empty()
    if (params.step == 'mapping' && !ref_bwa) {
        error "No BWA index configured. genome=${params.genome}, igenomes_base=${params.igenomes_base}, genome_bwa=${getGenomeAttribute('bwa')}"
    }
    ch_bwa = ref_bwa ? channel.fromPath(ref_bwa, checkIfExists: true).map { path -> [[id: 'bwa'], path] }.collect() : channel.empty()
    ch_chr_bed = ref_chr_bed ? channel.fromPath(ref_chr_bed).collect() : channel.empty()
    ch_centromere = ref_centromere ? channel.fromPath(ref_centromere).collect() : channel.value([])
    ch_medicc_arms = ref_medicc_arms ? channel.fromPath(ref_medicc_arms).collect() : channel.empty()
    ch_medicc_genes = ref_medicc_genes ? channel.fromPath(ref_medicc_genes).collect() : channel.empty()
    ch_gc_wig = ref_gc_wig ? channel.fromPath(ref_gc_wig).collect() : channel.empty()
    ch_map_wig = ref_map_wig_file ? channel.fromPath(ref_map_wig_file).collect() : channel.empty()
    ch_normal_wig = params.normal ? channel.fromPath(params.normal).collect() : channel.value([])
    //
    // SUBWORKFLOW: Run initialisation tasks
    //

    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.seq_center,
        params.seq_platform,
        params.step,
        params.library,
        ref_fasta
    )

    LP_WGS (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_fasta,
        ch_fasta_fai,
        ch_chr_arm_boundaries,
        ch_bwa,
        ch_chr_bed,
        ch_centromere,
        ch_medicc_arms,
        ch_medicc_genes,
        ch_gc_wig,
        ch_map_wig,
        ch_normal_wig,
        params.tools,
        params.genome,
        ref_qdnaseq_genome,
        ref_qdnaseq_package,
        ref_ichor_genome_build,
        ref_ichor_genome_style,
        params.step,
        params.tech,
        params.sort,
        params.fastp_adapter_fasta,
        params.filter_bam,
        params.filter_bam_min,
        params.filter_bam_max,
        params.call_gc,
        params.bin,
        params.ploidy,
        params.ascat_pcf_gamma,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        ch_bin_dir
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        LP_WGS.out
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
