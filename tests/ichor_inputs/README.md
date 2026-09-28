These routing tests run without containers or sequencing data, from the repository root:

```bash
nextflow -C /dev/null run tests/ichor_inputs/main.nf --scenario matched
nextflow -C /dev/null run tests/ichor_inputs/main.nf --scenario none
nextflow -C /dev/null run tests/ichor_inputs/main.nf --scenario duplicate
```

The first two must exit successfully with `PASS`. The duplicate scenario must
fail with `Multiple control WIGs for patient 'p1'`.

The samplesheet parser test can be run with Groovy and `lib` on the classpath:

```bash
groovy -cp lib tests/unit/samplesheet_status.groovy
```
