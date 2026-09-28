nextflow.enable.dsl = 2
include { ICHOR_INPUTS } from '../../subworkflows/local/ichor_inputs/main'

params.scenario = 'matched'
workflow {
    def records = [
        [[id: 'p1_t1', patient: 'p1', status: 1], 't1.wig'],
        [[id: 'p2_t2', patient: 'p2', status: 1], 't2.wig'],
        [[id: 'p1_t3', patient: 'p1', status: 1], 't3.wig']
    ]
    if (params.scenario in ['matched', 'duplicate']) {
        records << [[id: 'p1_n1', patient: 'p1', status: 0], 'n1.wig']
    }
    if (params.scenario == 'duplicate') {
        records << [[id: 'p1_n2', patient: 'p1', status: 0], 'n2.wig']
    }
    ICHOR_INPUTS(channel.fromList(records))
    ICHOR_INPUTS.out.inputs.toList().view { results ->
        assert results.size() == 3
        results.each { meta, tumour, normal ->
            assert meta.status == 1
            def expected = meta.patient == 'p1' && params.scenario == 'matched' ? 'n1.wig' : []
            assert normal == expected: "${meta.id}: ${normal} != ${expected}"
        }
        "PASS ${params.scenario}: all tumours retained and correctly paired"
    }
}
