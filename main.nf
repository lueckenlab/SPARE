nextflow.enable.dsl=2

// Nextflow requires a `main.nf` and a `nextflow.config` file
// to be present in the repository.
// Actual workflows are in the `src/workflows/` directory.
workflow {
    exit 0, "This is a dummy pipeline."
}