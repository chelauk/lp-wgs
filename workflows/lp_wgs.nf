/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { softwareVersionsToYAML      } from 'plugin/nf-core-utils'
include { MAPPING_QC                  } from '../subworkflows/local/mapping_qc/main'
include { CALLING_PREP                } from '../subworkflows/local/calling_prep/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
//include { REPORTING_MULTIQC           } from '../subworkflows/local/reporting_multiqc/main'
include { ICHORCNA_RUN                } from '../modules/local/ichorcna/run/main'
include { ICHORCNA_VERSIONS           } from '../modules/local/ichorcna/versions/main'
include { ACE                         } from '../modules/local/ace/main'
include { RUN_QDNASEQ                 } from '../modules/local/prep_ascat/main'
include { RUN_ASCAT                   } from '../modules/local/ascat_lp/main'
include { RUN_BAYES                   } from '../modules/local/bayes_cn/main'
include { PREP_MEDICC2                } from '../modules/local/prep_medicc2/main'
include { PREP_MEDICC2_ICHOR          } from '../modules/local/prep_medicc2_ichor/main'
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
    qdnaseq_package
    _ichor_genome_build
    _ichor_genome_style
    step
    tech
    sort
    fastp_adapter_fasta
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
    bin_dir

    main:
    selected_tools = tools.tokenize(',').collect { t -> t.trim() }.findAll { t -> t }    
    medicc_source = params.medicc_source ?: 'ace'

    if (qdnaseq_genome?.startsWith('mm')) {
        unsupported_tools = selected_tools.intersect(['medicc'])
        if (unsupported_tools) {
            exit 1, "Genome '${genome}' is configured as mouse (${qdnaseq_genome}), but these tools are still human-specific in this pipeline: ${unsupported_tools.join(', ')}."
        }
    }

    // define filter status
    filter_status = filter_bam ? "filter_${filter_bam_min}_${filter_bam_max}" : "filter_none"

    // To gather QC reports and software versions for reporting
    reports  = channel.empty()
    versions = channel.empty()

    // bin_dir for Rscripts
    bin_dir = channel.fromPath("$projectDir/bin").collect()
    ch_mapped_bam = channel.empty()

    if (step == 'mapping') {
        MAPPING_QC(
            ch_input_sample,
            bwa,
            fasta,
            fasta_fai,
            dict,
            chr_bed,
            sort,
            fastp_adapter_fasta,
            filter_bam,
            filter_bam_min,
            filter_bam_max,
            filter_status
        )
        ch_mapped_bam = MAPPING_QC.out.bam
    }

    CALLING_PREP(
        ch_input_sample,
        ch_mapped_bam,
        fasta,
        //fasta_fai,
        gc_wig,
        step,
        tech,
        filter_bam,
        filter_bam_min,
        filter_bam_max,
        call_gc,
    )
    ch_analysis_input = CALLING_PREP.out.analysis_input
    ch_gc_wig = CALLING_PREP.out.gc_wig
    versions = versions.mix(CALLING_PREP.out.versions)

    // run ichorcna
    if (selected_tools.contains('ichor')) {
        ICHORCNA_RUN(
            CALLING_PREP.out.readcounter_wig,
            ch_gc_wig,
            map_wig,
            normal_wig,
            [],
            centromere,
            [],
            []
        )
    }

    ICHORCNA_VERSIONS()

    // run QDNAseq once for ASCAT and/or ACE

    if (selected_tools.intersect(['ascat', 'ace'])) {
        RUN_QDNASEQ(ch_analysis_input, bin_size, qdnaseq_genome, qdnaseq_package)
    }

    if (selected_tools.contains('ascat')) {
        RUN_ASCAT(RUN_QDNASEQ.out.for_ascat, ploidy, chr_arm_boundaries, qdnaseq_genome, ascat_pcf_gamma)
    }

    // run ACE
    if (selected_tools.contains('ace')) {
        ACE(RUN_QDNASEQ.out.for_ace, filter_status, qdnaseq_genome, ploidy, bin_size)
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

    // run bayes_cna
    if (selected_tools.contains('bayes_cna')) {
        RUN_BAYES(ch_analysis_input, bin_size, qdnaseq_genome, bin_dir)
    }

    //run prep_medicc
    if (selected_tools.contains('medicc')) {
        if (medicc_source == 'ace') {
            if (!selected_tools.contains('ace')) {
                exit 1, "The 'medicc' workflow with medicc_source='ace' requires 'ace' so that ploidy-grouped inputs can be prepared."
            }
            PREP_MEDICC2(prep_medicc2_input, bin_dir, bin_size)
            versions = versions.mix(PREP_MEDICC2.out.versions)
            ch_medicc_input = PREP_MEDICC2.out.for_medicc
        } else if (medicc_source == 'ichor') {
            if (!selected_tools.contains('ichor')) {
                exit 1, "The 'medicc' workflow with medicc_source='ichor' requires 'ichor'."
            }
            ICHORCNA_RUN.out.cna_seg
                .map { meta, seg -> [meta.patient, meta.sample, meta.id, seg] }
                .groupTuple()
                .filter { tuple -> tuple[1].size() > 1 }
                .set { prep_medicc2_ichor_input }
            PREP_MEDICC2_ICHOR(prep_medicc2_ichor_input,bin_dir)
            versions = versions.mix(PREP_MEDICC2_ICHOR.out.versions)
            ch_medicc_input = PREP_MEDICC2_ICHOR.out.for_medicc
        } else {
            exit 1, "Unsupported medicc_source '${medicc_source}'. Supported values: ace, ichor."
        }

        // run medicc2
        MEDICC2(ch_medicc_input, medicc_arms, medicc_genes)
        versions = versions.mix(MEDICC2.out.versions)
    }

//    REPORTING_MULTIQC(
//        versions,
//        reports,
//        outdir,
//        multiqc_config,
//        multiqc_logo,
//        multiqc_methods_description
//    )
//
//    emit:
//    multiqc_report = REPORTING_MULTIQC.out // channel: /path/to/multiqc_report.html
//    versions
      def ch_multiqc_files = channel.empty()

      ch_multiqc_files = ch_multiqc_files
            .mix(MAPPING_QC.out.multiqc_files)
          //.mix(PROCESS_A.out.metrics.map { meta, file -> file })
          //.mix(PROCESS_B.out.report.map  { meta, file -> file })
      
      def ch_collated_versions = softwareVersionsToYAML(
          softwareVersions: channel.topic('versions'),
          nextflowVersion: workflow.nextflow.version,
      ).collectFile(
          storeDir: "${params.outdir}/pipeline_info",
          name: 'lp_wgs_software_mqc_versions.yml',
          sort: true,
          newLine: true,
      )
      
      ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

      MULTIQC(
             ch_multiqc_files
                 .flatten()
                 .collect()
                 .map { files ->
                     [
                         [id: 'lp_wgs'],
                         files,
                         multiqc_config
                             ? file(multiqc_config, checkIfExists: true)
                             : file(
                                 "${projectDir}/assets/multiqc_config.yml",
                                 checkIfExists: true
                             ),
                         [],
                         [],
                         [],
                     ]
                 }
         ) 
}
