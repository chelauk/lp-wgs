workflow SAMPLESHEET_TO_CHANNEL {
    take:
    ch_from_samplesheet
    seq_center
    seq_platform
    step

    main:
    
    ch_from_samplesheet
    /*
      .map { meta, fastq_1, fastq_2, bam, bai, rds ->
            if (fastq_1) {
                def flowcell = flowcellInfoFromFastq(fastq_1)
                def ID = "${flowcell[0]}:${flowcell[1]}:${flowcell[3]}"
                def PU = "${flowcell[2]}:${meta.patient}_${meta.sample}"
                meta += [
                    id         : "${meta.patient}_${meta.sample}${meta.lane ? "_${meta.lane}" : ""}",
                    read_group : "\"@RG\\tID:${ID}\\tCN:${seq_center}\\tPU:${PU}\\tSM:${meta.patient}_${meta.sample}\\tLB:${params.library}\\tDS:${params.fasta}\\tPL:${seq_platform}\"",
                    fastq      : true
                ]
                return [ meta, [ fastq_1, fastq_2 ] ]
            } else if (bam) {
                meta += [ bam: true, id: "${meta.patient}_${meta.sample}" ]
                if (rds) meta += [ rds: true ]
                meta.remove('lane')
				return [meta, *(rds ? [bam, bai, rds] : [bam, bai])] // spread operator.. need to read :-)
            }
        }
        .set { ch_samplesheet }
    */

      .map { t ->
      def (meta, fastq_1, fastq_2, bam, bai, rds) = t

      // normalise placeholders
      if (fastq_1 instanceof List && fastq_1.isEmpty()) fastq_1 = null
      if (fastq_2 instanceof List && fastq_2.isEmpty()) fastq_2 = null
      if (bam     instanceof List && bam.isEmpty())     bam     = null
      if (bai     instanceof List && bai.isEmpty())     bai     = null
      if (rds     instanceof List && rds.isEmpty())     rds     = null

      if (step == 'mapping') {
          if (!fastq_1 || !fastq_2)
              error "step=mapping requires fastq_1 and fastq_2 for ${meta.patient}_${meta.sample}"

                def flowcell = flowcellInfoFromFastq(fastq_1)
                def ID = "${flowcell[0]}:${flowcell[1]}:${flowcell[3]}"
                def PU = "${flowcell[2]}:${meta.patient}_${meta.sample}"
                meta = meta + [
                    id         : "${meta.patient}_${meta.sample}${meta.lane ? "_${meta.lane}" : ""}",
                    read_group : "\"@RG\\tID:${ID}\\tCN:${seq_center}\\tPU:${PU}\\tSM:${meta.patient}_${meta.sample}\\tLB:${params.library}\\tDS:${params.fasta}\\tPL:${seq_platform}\"",
                    data_type      : 'fastq'
                ]
          return tuple(meta, fastq_1, fastq_2)
      }

      else if (step == 'calling') {
          if (!bam)
              error "step=calling requires bam for ${meta.patient}_${meta.sample}"

          meta = meta + [
                 id       : "${meta.patient}_${meta.sample}",
                 data_type: 'bam'
               ]  
          return tuple(meta, bam, bai)
      }

      else if (step == 'tbd') {
          if (!rds)
              error "step=tbd requires rds for ${meta.patient}_${meta.sample}"

          meta = meta + [ data_type: 'rds' ]
          return tuple(meta, rds)
      }

      else {
          error "Unknown step '${step}'"
      }
    }
//    .view{ "samplesheet: $it" }
    .set{ ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def flowcellInfoFromFastq(file) {
    if (workflow.stubRun) {
        return [ "genericid", "genericnumber", "genericflowcell", "genericlane" ]
    }
    def line
    file.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    def fields = line.substring(1).split(':')
    return fields.size() >= 7 ?
        [ fields[0], fields[1], fields[2], fields[3] ] :
        [ "genericid", "genericnumber", fields[0], "genericlane" ]
}

