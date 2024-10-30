workflow SAMPLESHEET_TO_CHANNEL {
    take:
    ch_from_samplesheet             //
    seq_center                      //
    seq_platform                    //
    step                            //

    main:
    ch_from_samplesheet
        .map {
            meta, fastq_1, fastq_2, bam, bai, rds ->
            if ( fastq_1 ) {
                CN       = seq_center
                flowcell = flowcellInfoFromFastq(fastq_1)
                def ID = "${flowcell[0]}:${flowcell[1]}:${flowcell[3]}"
                def PU = "${flowcell[2]}:${meta.patient}_${meta.sample}"
                def read_group  = "\"@RG\\tID:${ID}\\tCN:${CN}\\tPU:${PU}\\tSM:${meta.patient}_${meta.sample}\\tLB:${params.library}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
                meta = meta + [ id: meta.patient + "_" + meta.sample + "_" + meta.lane ]
                meta = meta + [ rg: read_group ]
                return [ meta + [ fastq:true ], [ fastq_1, fastq_2 ] ]
                }
            else if ( bam ) {
			    if (rds) {
                    meta = meta + [ bam:true ]
				    meta = meta + [ rds:true ]
                    meta.remove('lane')
                    meta = meta + [ id: meta.patient + "_" + meta.sample ]
                    return [ meta, [ bam, bai, rds ]]
				} else {
                    meta = meta + [ bam:true ]
                    meta.remove('lane')
                    meta = meta + [ id: meta.patient + "_" + meta.sample ]
                    return [ meta, [ bam, bai ]]
					}
            }
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse first line of a FASTQ file, return the flowcell id and lane number.

def flowcellInfoFromFastq(file) {
    if (workflow.stubRun) {
        // Return default value when in -stub mode
        return "DEFAULT_FLOWCELL_LANE"
    }
    def line
    file.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String instrid = "genericid"
    String runid   = "genericnumber"
    String fcid    = "genericflowcell"
    String laneid  = "genericlane"

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        instrid = fields[0]
        runid   = fields[1]
        fcid    = fields[2]
        laneid  = fields[3]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return [ instrid, runid, fcid, laneid ]
}

/* Example usage with FASTQ header:
 instrument_id = "A02059"
run_number = "75"
flowcell_id = "H7GTYDRX5"
lane = "1"
sample_id = "Sample1"  # You can change this to a real sample name if known
library = "Library1"   # Descriptive library name
sequencing_center = "BI"
*/

