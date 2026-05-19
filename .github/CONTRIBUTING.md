# OceanOmics-OceanGenomes-Draft-Genomes: Contributing Guidelines

Hi there! Thanks for taking an interest in improving this pipeline.

Please use the issue and pull request templates where they help, but feel free to open a more general issue for discussion, planning, or support.

## Contribution workflow

1. Check whether there is already an issue or branch covering the change.
2. Make the necessary changes following the pipeline conventions below.
3. Run the relevant local tests.
4. Open a pull request against the development branch used by this repository.
5. Update `docs/usage.md`, `docs/output.md`, `README.md`, `CHANGELOG.md`, and `CITATIONS.md` where relevant.

## Tests

You can test changes locally by running the pipeline test profile:

```bash
nextflow run . -profile test,docker --outdir <OUTDIR>
```

For process selector warnings and other debug information, run:

```bash
nextflow run . -profile debug,test,docker --outdir <OUTDIR>
```

Module tests can be run with `nf-test` where test definitions are available:

```bash
nf-test test --profile debug,test,docker --verbose
```

## Pipeline contribution conventions

To make the code and processing logic easier to understand, keep changes close to the existing Nextflow structure and naming conventions.

### Adding a new step

1. Define the input channel into the new process from the expected upstream channel.
2. Write the process block.
3. Define the output channel if needed.
4. Add any new parameters to `nextflow.config` with sensible defaults.
5. Add new parameters to `nextflow_schema.json` with help text.
6. Add sanity checks and validation for relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. Update workflow tests where applicable.
9. Update `assets/multiqc_config.yml` if outputs should appear in MultiQC.
10. Document new output files in `docs/output.md`.

### Default values

Parameters should be defined with default values within the `params` scope in `nextflow.config`.

### Process resources

Sensible defaults for process resource requirements (CPUs, memory, and time) should be defined in `conf/base.config`, usually with `withLabel:` selectors that can be shared across multiple processes.

The process resources can be passed to tools dynamically with the `${task.cpus}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use naming schemes that make data flow easy to follow:

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Nextflow version

If you use a new core Nextflow feature, update the minimum required Nextflow version in `nextflow.config`.

### Images and figures

Keep overview images and documentation figures clear, lightweight, and specific to this pipeline.
