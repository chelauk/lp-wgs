/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MAPPING_QC                  } from '../subworkflows/local/mapping_qc/main'
include { CALLING_PREP                } from '../subworkflows/local/calling_prep/main'
include { REPORTING_MULTIQC           } from '../subworkflows/local/reporting_multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { ICHORCNA_RUN                } from '../modules/nf-core/ichorcna/run/main'
include { ACE                         } from '../modules/local/ace/main'
include { PREP_ASCAT                  } from '../modules/local/prep_ascat/main'
include { RUN_ASCAT                   } from '../modules/local/ascat_lp/main'
include { PREP_MEDICC2                } from '../modules/local/prep_medicc2/main'
include { MEDICC2                     } from '../modules/local/medicc2/main'


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
    normal_wig
    tools
    genome
    qdnaseq_genome
    step
    tech
    filter_bam
    filter_bam_min
    filter_bam_max
    call_gc
    bin_size
    ploidy
    ascat_pcf_gamma
    outdir
    multiqc_config
    multiqc_logo
    multiqc_methods_description

    main:
    selected_tools = tools.tokenize(',').collect { it.trim() }.findAll { it }

    if (qdnaseq_genome?.startsWith('mm')) {
        unsupported_tools = selected_tools.intersect(['medicc'])
        if (unsupported_tools) {
            exit 1, "Genome '${genome}' is configured as mouse (${qdnaseq_genome}), but these tools are still human-specific in this pipeline: ${unsupported_tools.join(', ')}."
        }
    }

    // define filter status
    filter_status = filter_bam ? "filter_${filter_bam_min}_${filter_bam_max}" : "filter_none"

    // To gather QC reports and software versions for reporting
    reports  = Channel.empty()
    versions = Channel.empty()

    // bin_dir for Rscripts
    bin_dir = Channel.fromPath("$projectDir/bin").collect()

    if (step == 'mapping') {
        MAPPING_QC(
            ch_input_sample,
            bwa,
            fasta,
            dict,
            chr_bed,
            filter_bam,
            filter_bam_min,
            filter_bam_max,
            filter_status
        )
        ch_mapped_bam = MAPPING_QC.out.bam
        versions = versions.mix(MAPPING_QC.out.versions)
        reports  = reports.mix(MAPPING_QC.out.reports)
    }

    CALLING_PREP(
        ch_input_sample,
        ch_mapped_bam,
        fasta,
        gc_wig,
        step,
        tech,
        filter_bam,
        filter_bam_min,
        filter_bam_max,
        call_gc,
        bin_size
    )
    ch_analysis_input = CALLING_PREP.out.analysis_input
    ch_gc_wig = CALLING_PREP.out.gc_wig
    versions = versions.mix(CALLING_PREP.out.versions)

    // run ichorcna
    if (selected_tools.contains('ichor')) {
        ICHORCNA_RUN(
            CALLING_PREP.out.readcounter_wig,
            normal_wig,
            ch_gc_wig,
            map_wig,
            centromere,
            filter_status
        )
        versions= versions.mix(ICHORCNA_RUN.out.versions)
    }

    // run PREP_ASCAT

    if (selected_tools.contains('ascat')) {
        PREP_ASCAT(ch_analysis_input, bin_size)
        RUN_ASCAT(PREP_ASCAT.out.for_ascat, ploidy, chr_arm_boundaries, ascat_pcf_gamma)
    }

    // run ACE
    if (selected_tools.contains('ace')) {
        ACE(ch_analysis_input, filter_status)
        versions = versions.mix(ACE.out.versions)
        ACE.out.ace
            .map { meta, ace ->
                // If meta.predicted_ploidy is null, set it to 2
                meta.predicted_ploidy = meta.predicted_ploidy ?: 2
                [meta.patient, meta.sample, meta.id, meta.predicted_ploidy, ace]
            }
            .groupTuple()
            .filter { tuple -> tuple[1].size() > 1 }
            .set { prep_medicc2_input }
    }

    //run prep_medicc
    if (selected_tools.contains('medicc')) {
        if (!selected_tools.contains('ace')) {
            exit 1, "The 'medicc' workflow currently requires 'ace' so that ploidy-grouped inputs can be prepared."
        }
        PREP_MEDICC2(prep_medicc2_input, bin_dir)
        versions = versions.mix(PREP_MEDICC2.out.versions)

        // run medicc2
        MEDICC2(PREP_MEDICC2.out.for_medicc, medicc_arms, medicc_genes)
    }

    REPORTING_MULTIQC(
        versions,
        reports,
        outdir,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description
    )

    emit:
    multiqc_report = REPORTING_MULTIQC.out.report // channel: /path/to/multiqc_report.html
    versions
}
