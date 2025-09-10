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

params.dict                  = getGenomeAttribute('dict')
params.fasta                 = getGenomeAttribute('fasta')
params.fasta_fai             = getGenomeAttribute('fasta_fai')
params.bwa                   = getGenomeAttribute('bwa')
params.centromere            = getGenomeAttribute('centromere')
params.map_wig               = getGenomeAttribute('map_wig')
params.chr_bed               = getGenomeAttribute('chr_bed')
params.medicc_arms           = getGenomeAttribute('medicc_arms')
params.medicc_genes          = getGenomeAttribute('medicc_genes')
params.chr_arm_boundaries    = getGenomeAttribute('chr_arm_boundaries')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_lp_wgs_pipeline'
include { LP_WGS                           } from './workflows/lp_wgs'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_lp_wgs_pipeline'


// Initialise fasta file with meta map:
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it] }.collect(): Channel.empty


//
// gather prebuilt indices
//
dict                   = params.dict               ? Channel.fromPath(params.dict).collect()               : Channel.empty()
fasta_fai              = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it] }.collect()          : Channel.empty()
chr_arm_boundaries     = params.chr_arm_boundaries ? Channel.fromPath(params.chr_arm_boundaries).collect() : Channel.empty()
bwa                    = params.bwa                ? Channel.fromPath(params.bwa).map{ it -> [ [id:'bwa'], it] }.collect()                : Channel.empty()
chr_bed                = params.chr_bed            ? Channel.fromPath(params.chr_bed).collect()            : Channel.empty()
centromere             = params.centromere         ? Channel.fromPath(params.centromere).collect()         : Channel.empty()
medicc_arms            = params.medicc_arms        ? Channel.fromPath(params.medicc_arms).collect()        : Channel.empty()
medicc_genes           = params.medicc_genes       ? Channel.fromPath(params.medicc_genes).collect()       : Channel.empty()
gc_wig                 = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_${params.bin}kb.wig").collect()   : Channel.empty()
map_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_${params.bin}kb.wig").collect()   : Channel.empty()
pon_rds                = params.map_wig            ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_hg38_${params.bin}kb_median_normAutosome_median.rds").collect() : Channel.value([]) // optional Channel.empty()
normal_wig             = params.normal_wig         ? Channel.fromPath(params.normal_wig).collect()         : Channel.value([]) // empty value channel necessary

//
// WORKFLOW: Run main lp-wgs analysis pipeline
//
workflow {
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
        params.input
    )

    LP_WGS (
        PIPELINE_INITIALISATION.out.samplesheet,
        fasta,
        dict,
        fasta_fai,
        chr_arm_boundaries,
        bwa,
        chr_bed,
        centromere,
        medicc_arms,
        medicc_genes,
        gc_wig,
        map_wig,
        pon_rds,
        normal_wig
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
		if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
