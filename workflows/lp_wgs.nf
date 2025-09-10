/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_lp_wgs_pipeline'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.bwa,
    params.fasta,
    params.fasta_fai,
    params.centromere,
    params.chr_arm_boundaries,
    params.map_wig,
    params.dict,
    params.medicc_arms,
    params.medicc_genes
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_COMPLETION              } from '../subworkflows/local/utils_nfcore_lp_wgs_pipeline'
include { PIPELINE_INITIALISATION          } from '../subworkflows/local/utils_nfcore_lp_wgs_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { QC_TRIM                     } from '../subworkflows/local/qc_trim/main'
include { MERGE_LANES                 } from '../subworkflows/local/merge_lanes/main'
include { MOSDEPTH                    } from '../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { HMMCOPY_GCCOUNTER           } from '../modules/nf-core/hmmcopy/gccounter/main'
include { HMMCOPY_READCOUNTER         } from '../modules/nf-core/hmmcopy/readcounter/main'
include { ICHORCNA_RUN                } from '../modules/nf-core/ichorcna/run/main'
include { SAMTOOLS_VIEW               } from '../modules/nf-core/samtools/view/'
include { SAMTOOLS_NVIEW              } from '../modules/nf-core/samtools/view_nanopore/main'
include { ACE                         } from '../modules/local/ace'
include { PREP_ASCAT                  } from '../modules/local/prep_ascat/main'
include { RUN_ASCAT                   } from '../modules/local/ascat_lp/main'
include { PREP_MEDICC2                } from '../modules/local/prep_medicc2/main'
include { MEDICC2                     } from '../modules/local/medicc2/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../modules/local/picard/collectalignmentmummarymetrics/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LP_WGS {

    take:
    ch_input_sample
    fasta
    dict
    fasta_fai
    chr_arm_boundaries
    bwa
    chr_bed
    centromere
    medicc_arms
    medicc_genes
    gc_wig
    map_wig
    pon_rds
    normal_wig

    main:
    // define filter status
    filter_status = params.filter_bam ? "filter_${params.filter_bam_min}_${params.filter_bam_max}" : "filter_none"

    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    reports          = Channel.empty()
    versions         = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // To gather mosdepth reports
    // mosdepth_reports = Channel.empty()

    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    //
    // FASTQC, FASTP and BWA
    //
    sort_bam = true
    // Gather index for mapping given the chosen aligner
    ch_map_index = bwa
    // bin_dir for Rscripts
    bin_dir = Channel.fromPath("$projectDir/bin").collect()

        // Mapping step
    println( params.step )
    if ( params.step == 'mapping' ) {
        // 1. QC and Trim with fastqc and fastp ( QC_TRIM subworkflow )
        QC_TRIM ( ch_input_sample, ch_map_index, sort_bam)
        versions = versions.mix(QC_TRIM.out.versions)
        reports  = reports.mix(QC_TRIM.out.reports)

        // 2. Map with BWA MEM
        BWA_MEM( QC_TRIM.out.reads, ch_map_index, sort_bam) // If aligner is bwa-mem
        versions = versions.mix(BWA_MEM.out.versions.first())

        // 3. call the merge subworkflow
        MERGE_LANES ( BWA_MEM.out.bam )
        versions = versions.mix(MERGE_LANES.out.versions.first())

        // 4. filter bams
        if ( !params.filter_bam ) {
            ch_bam_input = MERGE_LANES.out.bam
            } else {
                ch_filter_input = MERGE_LANES.out.bam
                SAMTOOLS_VIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
                ch_bam_input = SAMTOOLS_VIEW.out.bam
                versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
            }
        
        PICARD_COLLECTALIGNMENTSUMMARYMETRICS (
            ch_bam_input,
            fasta,
            dict,
            filter_status
            )
        versions = versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions.first())
        reports  = reports.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect{meta, report -> report})
        PICARD_COLLECTINSERTSIZEMETRICS (
            ch_bam_input,
            filter_status
            )
        versions = versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())
        reports  = reports.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.size_metrics.collect{meta, report -> report}) 
      
        MOSDEPTH(
            ch_bam_input,
            chr_bed,
            fasta,
            filter_status
        )
        //reports  = reports.mix(MOSDEPTH.out.reports.collect{meta, report -> report})
        versions = versions.mix(MOSDEPTH.out.versions.first())
    }
    
    // Calling step

    else if ( params.step == 'calling' ) {
        ch_bam_input = ch_input_sample
        // 1. run the illumina tech branch
        if ( params.tech == "illumina"){
            ch_bam_input = ch_input_sample
                            .map { meta, files ->
                                  [meta, files[0], files[1]] }
            if ( params.filter_bam ) {
               ch_filter_input = ch_bam_input
               SAMTOOLS_VIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
               ch_bam_input = SAMTOOLS_VIEW.out.bam
               versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
            }


            } else if (params.tech == "nanopore"){
                // 2. run the nanopore tech branch
                ch_bam_input = ch_input_sample
                                .map { meta, files ->
                                [meta, files[0], files[1]] }
                if ( params.filter_bam ) {
                    ch_filter_input = ch_bam_input
                        SAMTOOLS_NVIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
                        ch_bam_input = SAMTOOLS_NVIEW.out.bam
                                         .map{ meta, bam, bai -> [ meta, [ bam, bai]] }
                        versions = versions.mix(SAMTOOLS_NVIEW.out.versions.first())
                     }
            }
    }
    // run hmmcopygccounter
    if ( params.call_gc ) {
        HMMCOPY_GCCOUNTER(fasta, params.bin )
        gc_wig = HMMCOPY_GCCOUNTER.out.wig
        versions = versions.mix(HMMCOPY_GCCOUNTER.out.versions)
        } else {
            gc_wig = gc_wig
        }

    HMMCOPY_READCOUNTER( ch_bam_input )
    versions = versions.mix(HMMCOPY_READCOUNTER.out.versions)

    // run ichorcna
    if  ( params.step != 'ascat' ) {
        ICHORCNA_RUN(
            HMMCOPY_READCOUNTER.out.wig,
            normal_wig,
            gc_wig,
            map_wig,
            pon_rds,
            centromere,
            filter_status
        )
        versions= versions.mix(ICHORCNA_RUN.out.versions)
    }

    // run PREP_ASCAT

    if (params.step == 'ascat') {
        println("bin " + params.bin)
        PREP_ASCAT( ch_bam_input, params.bin )
        RUN_ASCAT( PREP_ASCAT.out.for_ascat, params.ploidy, chr_arm_boundaries )
    }

    // run ACE
    if (params.step != 'ascat') {
    ACE(ch_bam_input, filter_status)
    versions = versions.mix(ACE.out.versions)

    ACE.out.ace
        .map{ meta, ace ->
        // If meta.predicted_ploidy is null, set it to 2
        meta.predicted_ploidy = meta.predicted_ploidy ?: 2
        [meta.patient, meta.sample, meta.id, meta.predicted_ploidy, ace]
            }
        .groupTuple()
        .filter { tuple -> tuple[1].size() > 1 }
        .set{ prep_medicc2_input }
    }
    
    if ( params.step == 'ascat' ) {
    prep_medicc2_input = ch_bam_input
    }
    
    //run prep_medicc
    PREP_MEDICC2(prep_medicc2_input, bin_dir)
       versions = versions.mix(PREP_MEDICC2.out.versions)

    // run medicc2
    MEDICC2(PREP_MEDICC2.out.for_medicc, medicc_arms, medicc_genes)

    //
    // Collate and save software versions
    //
    version_yaml = Channel.empty()
    version_yaml = softwareVersionsToYAML(versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_sarek_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(version_yaml)
    ch_multiqc_files                      = ch_multiqc_files.mix(reports)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))



    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    emit:
        multiqc_report // channel: /path/to/multiqc_report.html
        versions
}
