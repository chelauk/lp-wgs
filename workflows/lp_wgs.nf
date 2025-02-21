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

//
// gather prebuilt indices
//
    dict                   = params.dict               ? Channel.fromPath(params.dict).collect()               : Channel.empty()
    fasta                  = params.fasta              ? Channel.fromPath(params.fasta).collect()              : Channel.empty()
	fasta_fai              = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()          : Channel.empty()
    chr_arm_boundaries     = params.chr_arm_boundaries ? Channel.fromPath(params.chr_arm_boundaries).collect() : Channel.empty()
    bwa                    = params.bwa                ? Channel.fromPath(params.bwa).collect()                : Channel.empty()
    chr_bed                = params.chr_bed            ? Channel.fromPath(params.chr_bed).collect()            : Channel.empty()
    centromere             = params.centromere         ? Channel.fromPath(params.centromere).collect()         : Channel.empty()
    medicc_arms            = params.medicc_arms        ? Channel.fromPath(params.medicc_arms).collect()        : Channel.empty()
    medicc_genes           = params.medicc_genes       ? Channel.fromPath(params.medicc_genes).collect()       : Channel.empty()
    gc_wig                 = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_${params.bin}kb.wig").collect()   : Channel.empty()
	map_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_${params.bin}kb.wig").collect()   : Channel.empty()
	pon_rds                = params.map_wig            ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_hg38_${params.bin}kb_median_normAutosome_median.rds").collect() : Channel.value([]) // optional Channel.empty()
    normal_wig             = params.normal_wig         ? Channel.fromPath(params.normal_wig).collect()         : Channel.value([]) // empty value channel necessary

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LP_WGS {

    take:
    ch_input_sample

    main:
    // define filter status
    if ( params.step != 'ascat' ) {
	filter_status = params.filter_bam == null ? "filter_none" : "filter_" + params.filter_bam_min + "_" + params.filter_bam_max
	}

    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    reports          = Channel.empty()
    versions         = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // To gather mosdepth reports
    mosdepth_reports = Channel.empty()

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

    if ( params.step == 'fastq' ) {
        fastq_input = ch_input_sample
        QC_TRIM ( fastq_input, ch_map_index, sort_bam)
        versions = versions.mix(QC_TRIM.out.versions)
        reports  = reports.mix(QC_TRIM.out.reports)
        QC_TRIM.out.reads.view{ "Trimmed reads: $it" }
        BWA_MEM( QC_TRIM.out.reads, ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort_bam) // If aligner is bwa-mem
        versions = versions.mix(BWA_MEM.out.versions.first())
        MERGE_LANES ( BWA_MEM.out.bam )
        versions = versions.mix(MERGE_LANES.out.versions.first())
        if ( params.filter_bam == null ) {
            ch_bam_input = MERGE_LANES.out.bam
            } else if ( params.filter_bam != null  ) {
                ch_filter_input = MERGE_LANES.out.bam
                SAMTOOLS_VIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
                ch_bam_input = SAMTOOLS_VIEW.out.bam
                versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
            }
        }
    if ( params.step == 'bam' || params.step == 'ascat' ) {
        ch_bam_input = ch_input_sample
	}
    if ( params.step == 'bam'  &&  params.filter_bam != null && params.tech != "nanopore"){
        ch_filter_input = ch_input_sample
                            .map { meta, files ->
							    [meta, files[0], files[1]] }
        SAMTOOLS_VIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
        ch_bam_input = SAMTOOLS_VIEW.out.bam
        versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
    }
    if ( params.step == 'bam'  &&  params.filter_bam != null && params.tech == "nanopore"){
        ch_filter_input = ch_input_sample
                            .map { meta, files ->
							    [meta, files[0], files[1]] }
        SAMTOOLS_NVIEW ( ch_filter_input, params.filter_bam_min, params.filter_bam_max )
        ch_bam_input = SAMTOOLS_NVIEW.out.bam
		                            .map{ meta, bam, bai -> [ meta, [ bam, bai]] }
		versions = versions.mix(SAMTOOLS_NVIEW.out.versions.first())
    }
    if (( params.step != 'ascat' ) && ( params.tech == 'illumina' )) {
		PICARD_COLLECTALIGNMENTSUMMARYMETRICS ( ch_bam_input , fasta, dict,filter_status)
        versions = versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions.first())
        reports  = reports.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect{meta, report -> report})
		PICARD_COLLECTINSERTSIZEMETRICS ( ch_bam_input ,filter_status)
        versions = versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())
        reports  = reports.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.size_metrics.collect{meta, report -> report})

		MOSDEPTH(
            ch_bam_input,
            chr_bed,
            fasta.map{ it -> [[id:it[0].baseName], it] },
            filter_status
        )
        mosdepth_reports = mosdepth_reports.mix(MOSDEPTH.out.global_txt,
                                        MOSDEPTH.out.regions_txt)
        reports  = reports.mix(mosdepth_reports.collect{meta, report -> report})
        versions = versions.mix(MOSDEPTH.out.versions.first())
    }

    // run hmmcopygccounter
    if ( params.call_gc ) {
        HMMCOPY_GCCOUNTER(fasta.map{ it -> [[id:it[0].baseName], it] }, params.bin )
        gc_wig = HMMCOPY_GCCOUNTER.out.wig
        versions = versions.mix(HMMCOPY_GCCOUNTER.out.versions)
        } else {
            gc_wig = gc_wig
        }

    // run hmmcopyreadcounter
    if ( params.step != 'ascat' && params.step != 'bam') {
        HMMCOPY_READCOUNTER( ch_bam_input )
        versions = versions.mix(HMMCOPY_READCOUNTER.out.versions)
    } else if ( params.step == 'bam') {
	    ch_bam_input = ch_bam_input
		                .map { meta, files -> [ meta, files[0], files[1]]}
        HMMCOPY_READCOUNTER( ch_bam_input )
        versions = versions.mix(HMMCOPY_READCOUNTER.out.versions)
    }

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
	ch_bam_input
		.map{ meta , files ->
		 meta = meta + [ id: meta.patient + "_" + meta.sample ]
		 // If meta.predicted_ploidy is null, set it to 2
		 meta.predicted_ploidy = meta.predicted_ploidy ?: 2
		 [meta.patient, meta.sample, meta.id, meta.predicted_ploidy, files]
		 }
		.groupTuple()
		.map { patient, sample, id, ploidy, files ->
		  files = files.flatten()
		  [patient, sample, id, ploidy, files ]
		  }
		.filter { tuple -> tuple[1].size() > 1 }
		.set{ prep_medicc2_input }
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
