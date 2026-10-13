// Keep each tumour and its selected normal in one tuple so asynchronous
// readCounter tasks cannot pair a tumour with another patient's control.
workflow ICHOR_INPUTS {
    take:
    wigs

    main:
    controls = wigs
        .filter { meta, _wig -> meta.status == 0 }
        .map { meta, wig -> [meta.patient, wig] }
        .toList()
        .map { records ->
            def normals = [:]
            records.each { patient, wig ->
                if (normals.containsKey(patient)) {
                    error "Multiple control WIGs for patient '${patient}'. Supply one control sample per patient."
                }
                normals[patient] = wig
            }
            normals
        }

    emit:
    wigs
        .filter { meta, _wig -> meta.status != 0 }
        .combine(controls)
        .map { meta, wig, normals ->
            tuple(meta, wig, normals.containsKey(meta.patient) ? normals[meta.patient] : [])
        }
}
