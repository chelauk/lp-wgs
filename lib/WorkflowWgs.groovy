//
// This file holds several functions specific to the workflow/wgs.nf in the lp-wgs pipeline
//

import groovy.text.SimpleTemplateEngine
import groovy.json.JsonSlurper

class WorkflowWgs {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static LinkedHashMap paramsSummaryMap(workflow, params, schema_filename='nextflow_schema.json') {
        def workflow_summary = [:]
        if (workflow.revision) {
            workflow_summary['revision'] = workflow.revision
        }
        workflow_summary['runName'] = workflow.runName
        if (workflow.containerEngine) {
            workflow_summary['containerEngine'] = workflow.containerEngine
        }
        if (workflow.container) {
            workflow_summary['container'] = workflow.container
        }
        workflow_summary['launchDir'] = workflow.launchDir
        workflow_summary['workDir'] = workflow.workDir
        workflow_summary['projectDir'] = workflow.projectDir
        workflow_summary['userName'] = workflow.userName
        workflow_summary['profile'] = workflow.profile
        workflow_summary['configFiles'] = workflow.configFiles.join(', ')

        def schema = new JsonSlurper().parse(new File("${workflow.projectDir}/${schema_filename}"))
        def groups = []
        schema.allOf?.each { ref ->
            def ref_path = ref['$ref']
            if (ref_path?.startsWith('#/$defs/')) {
                def group_name = ref_path.replace('#/$defs/', '')
                def group = schema['$defs']?.get(group_name)
                if (group) {
                    groups << group
                }
            }
        }
        if (schema.properties) {
            groups << [title: 'Other parameters', properties: schema.properties]
        }

        def params_summary = [:]
        groups.each { group ->
            def title = group.title ?: 'Parameters'
            def sub_params = new LinkedHashMap()
            group.properties?.keySet()?.sort()?.each { param ->
                if (params.containsKey(param)) {
                    def params_value = params.get(param)
                    def schema_value = group.properties[param].default
                    if (schema_value != null && params_value != schema_value) {
                        sub_params.put(param, params_value)
                    } else if (schema_value == null && params_value != '' && params_value != null && params_value != false) {
                        sub_params.put(param, params_value)
                    }
                }
            }
            params_summary.put(title, sub_params)
        }
        return ['Core Nextflow options': workflow_summary] << params_summary
    }

    public static List parseSamplesheet(String samplesheet) {
        def lines = new File(samplesheet).readLines().findAll { it.trim() }
        if (!lines) {
            return []
        }
        def header = splitCsvLine(lines[0])
        return lines.drop(1).collect { line ->
            def values = splitCsvLine(line)
            def row = [:]
            header.eachWithIndex { name, index ->
                row[name] = index < values.size() ? values[index] : ''
            }
            row
        }
    }

    private static List splitCsvLine(String line) {
        def values = []
        def current = new StringBuilder()
        boolean in_quotes = false
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i)
            if (c == '"') {
                in_quotes = !in_quotes
            } else if (c == ',' && !in_quotes) {
                values << current.toString()
                current = new StringBuilder()
            } else {
                current.append(c)
            }
        }
        values << current.toString()
        return values.collect { it.trim() }
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }}
