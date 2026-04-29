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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.dict                  = params.dict ?: getGenomeAttribute('dict')
params.fasta                 = params.fasta ?: getGenomeAttribute('fasta')
params.fasta_fai             = params.fasta_fai ?: getGenomeAttribute('fasta_fai')
params.bwa                   = params.bwa ?: getGenomeAttribute('bwa')
params.centromere            = params.centromere ?: getGenomeAttribute('centromere')
params.map_wig               = params.map_wig ?: getGenomeAttribute('map_wig')
params.map_wig_file          = params.map_wig_file ?: getGenomeAttribute('map_wig_file')
params.gc_wig                = params.gc_wig ?: getGenomeAttribute('gc_wig')
params.ichor_genome_build    = params.ichor_genome_build ?: getGenomeAttribute('ichor_genome_build')
params.ichor_genome_style    = params.ichor_genome_style ?: getGenomeAttribute('ichor_genome_style')
params.chr_bed               = params.chr_bed ?: getGenomeAttribute('chr_bed')
params.medicc_arms           = params.medicc_arms ?: getGenomeAttribute('medicc_arms')
params.medicc_genes          = params.medicc_genes ?: getGenomeAttribute('medicc_genes')
params.chr_arm_boundaries    = params.chr_arm_boundaries ?: getGenomeAttribute('chr_arm_boundaries')
params.qdnaseq_genome        = params.qdnaseq_genome ?: getGenomeAttribute('qdnaseq_genome')
params.qdnaseq_package       = params.qdnaseq_package ?: getGenomeAttribute('qdnaseq_package')
params.hmmcopy_chromosomes   = params.hmmcopy_chromosomes ?: getGenomeAttribute('hmmcopy_chromosomes')

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
    main:
    // Initialise genome resources close to the workflow entrypoint.
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta).map { path -> [[id: path.baseName], path] }.collect() : Channel.empty()

    ch_dict = params.dict ? Channel.fromPath(params.dict).collect() : Channel.empty()
    ch_fasta_fai = params.fasta_fai ? Channel.fromPath(params.fasta_fai).map { path -> [[id: 'fai'], path] }.collect() : Channel.empty()
    ch_chr_arm_boundaries = params.chr_arm_boundaries ? Channel.fromPath(params.chr_arm_boundaries).collect() : Channel.empty()
    ch_bwa = params.bwa ? Channel.fromPath(params.bwa).map { path -> [[id: 'bwa'], path] }.collect() : Channel.empty()
    ch_chr_bed = params.chr_bed ? Channel.fromPath(params.chr_bed).collect() : Channel.empty()
    ch_centromere = params.centromere ? Channel.fromPath(params.centromere).collect() : Channel.value([])
    ch_medicc_arms = params.medicc_arms ? Channel.fromPath(params.medicc_arms).collect() : Channel.empty()
    ch_medicc_genes = params.medicc_genes ? Channel.fromPath(params.medicc_genes).collect() : Channel.empty()
    ch_gc_wig = params.gc_wig ? Channel.fromPath(params.gc_wig).collect() : Channel.empty()
    ch_map_wig = params.map_wig_file ? Channel.fromPath(params.map_wig_file).collect() : Channel.empty()
    ch_normal_wig = params.normal ? Channel.fromPath(params.normal).collect() : Channel.value([])

    //
    // SUBWORKFLOW: Run initialisation tasks
    //

    PIPELINE_INITIALISATION(
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.seq_center,
        params.seq_platform,
        params.step
    )

    LP_WGS (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_fasta,
        ch_dict,
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
        params.qdnaseq_genome,
        params.step,
        params.tech,
        params.filter_bam,
        params.filter_bam_min,
        params.filter_bam_max,
        params.call_gc,
        params.bin,
        params.ploidy,
        params.outdir,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description
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
        LP_WGS.out.multiqc_report
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
