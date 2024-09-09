workflow SAMPLESHEET_TO_CHANNEL {
    take:
    ch_from_samplesheet             //
    seq_center                      //
    seq_platform                    //
    step                            //

   main:
    ch_from_samplesheet
        .map {
            meta, fastq_1, fastq_2, bam, bai ->
            if ( fastq_1 ) {
                meta = meta + [ id: meta.patient + "_" + meta.sample + "_" + meta.lane ] 
                return [ meta + [ fastq:true ], [ fastq_1, fastq_2 ] ]
                }
            else if ( bam ) {
                meta = meta + [ bam:true ] 
                meta.remove('lane')
                meta = meta + [ id: meta.patient + "_" + meta.sample ] 
                return [ meta, [ bam, bai ]]
            }
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
}
