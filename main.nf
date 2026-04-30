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

ref_dict                  = params.dict ?: getGenomeAttribute('dict')
ref_fasta                 = params.fasta ?: getGenomeAttribute('fasta')
ref_fasta_fai             = params.fasta_fai ?: getGenomeAttribute('fasta_fai')
ref_bwa                   = params.bwa ?: getGenomeAttribute('bwa')
ref_centromere            = params.centromere ?: getGenomeAttribute('centromere')
ref_map_wig               = params.map_wig ?: getGenomeAttribute('map_wig')
ref_map_wig_file          = params.map_wig_file ?: getGenomeAttribute('map_wig_file')
ref_gc_wig                = params.gc_wig ?: getGenomeAttribute('gc_wig')
ref_ichor_genome_build    = params.ichor_genome_build ?: getGenomeAttribute('ichor_genome_build')
ref_ichor_genome_style    = params.ichor_genome_style ?: getGenomeAttribute('ichor_genome_style')
ref_chr_bed               = params.chr_bed ?: getGenomeAttribute('chr_bed')
ref_medicc_arms           = params.medicc_arms ?: getGenomeAttribute('medicc_arms')
ref_medicc_genes          = params.medicc_genes ?: getGenomeAttribute('medicc_genes')
ref_chr_arm_boundaries    = params.chr_arm_boundaries ?: getGenomeAttribute('chr_arm_boundaries')
ref_qdnaseq_genome        = params.qdnaseq_genome ?: getGenomeAttribute('qdnaseq_genome')
ref_qdnaseq_package       = params.qdnaseq_package ?: getGenomeAttribute('qdnaseq_package')
ref_hmmcopy_chromosomes   = params.hmmcopy_chromosomes ?: getGenomeAttribute('hmmcopy_chromosomes')

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
    ch_fasta = ref_fasta ? Channel.fromPath(ref_fasta).map { path -> [[id: path.baseName], path] }.collect() : Channel.empty()

    ch_dict = ref_dict ? Channel.fromPath(ref_dict).collect() : Channel.empty()
    ch_fasta_fai = ref_fasta_fai ? Channel.fromPath(ref_fasta_fai).map { path -> [[id: 'fai'], path] }.collect() : Channel.empty()
    ch_chr_arm_boundaries = ref_chr_arm_boundaries ? Channel.fromPath(ref_chr_arm_boundaries).collect() : Channel.empty()
    if (params.step == 'mapping' && !ref_bwa) {
        error "No BWA index configured. genome=${params.genome}, igenomes_base=${params.igenomes_base}, genome_bwa=${getGenomeAttribute('bwa')}"
    }
    ch_bwa = ref_bwa ? Channel.fromPath(ref_bwa, checkIfExists: true).map { path -> [[id: 'bwa'], path] }.collect() : Channel.empty()
    ch_chr_bed = ref_chr_bed ? Channel.fromPath(ref_chr_bed).collect() : Channel.empty()
    ch_centromere = ref_centromere ? Channel.fromPath(ref_centromere).collect() : Channel.value([])
    ch_medicc_arms = ref_medicc_arms ? Channel.fromPath(ref_medicc_arms).collect() : Channel.empty()
    ch_medicc_genes = ref_medicc_genes ? Channel.fromPath(ref_medicc_genes).collect() : Channel.empty()
    ch_gc_wig = ref_gc_wig ? Channel.fromPath(ref_gc_wig).collect() : Channel.empty()
    ch_map_wig = ref_map_wig_file ? Channel.fromPath(ref_map_wig_file).collect() : Channel.empty()
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
        params.step,
        params.library,
        ref_fasta
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
        ref_qdnaseq_genome,
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
