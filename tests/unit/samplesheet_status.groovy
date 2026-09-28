def checkSheet(String content, Closure check) {
    def sheet = File.createTempFile('lp-wgs-status-', '.csv')
    try {
        sheet.text = content
        check(WorkflowWgs.parseSamplesheet(sheet.path))
    } finally {
        sheet.delete()
    }
}

checkSheet('patient,sample\np1,t1\n') { assert it[0].status == 1 }
checkSheet('patient,sample,status\np1,t1,1\np1,n1,0\np2,t2,\n') {
    assert it*.status == [1, 0, 1]
}
checkSheet('patient,sample,lane,status\np1,n1,1,0\np1,n1,2,0\n') {
    assert it*.status == [0, 0]
}
[
    'patient,sample,status\np1,t1,2\n',
    'patient,sample,status\np1,t1,tumour\n',
    'patient,sample,lane,status\np1,t1,1,0\np1,t1,2,1\n'
].each { content ->
    try {
        checkSheet(content) { }
        assert false: 'Expected invalid status to fail'
    } catch (IllegalArgumentException expected) {
        assert expected.message.contains('status')
    }
}
println 'Samplesheet status tests passed'
