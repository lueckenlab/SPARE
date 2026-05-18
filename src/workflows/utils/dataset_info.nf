// Loader for per-dataset config files. See dataset_info_example.yaml
// for the schema.
//
// Usage in a workflow:
//   include { expandDatasetInfo } from "${params.rootDir}/src/workflows/utils/dataset_info.nf"
//   input_ch | map { id, state ->
//     def info = readYaml(state.dataset_info)
//     [id, state + expandDatasetInfo(info, state.dataset_info)]
//   }
//
// Returns a flat map keyed to the argument names each Viash component
// already expects, so per-stage `fromState` blocks stay short.

// Resolve a file referenced inside dataset_info.yaml against the data
// directory that holds the YAML itself.
def resolveDataPath(String yamlPath, relPath) {
    if (relPath == null) return null
    return new File(new File(yamlPath).parentFile, relPath as String).absolutePath
}

def expandDatasetInfo(Map info, String yamlPath) {
    def files = info.files ?: [:]
    def keys = info.keys ?: [:]
    def meta = info.metadata ?: [:]
    def pp = info.preprocess ?: [:]
    def ev = info.evaluate ?: [:]

    return [
        // file paths (absolute, resolved against data/<id>/)
        raw_path:             resolveDataPath(yamlPath, files.raw),
        cleaned_path:         resolveDataPath(yamlPath, files.cleaned),
        processed_path:       resolveDataPath(yamlPath, files.processed),
        metadata_path:        resolveDataPath(yamlPath, files.metadata),
        representations_path: resolveDataPath(yamlPath, files.representations),
        benchmark_schema:     resolveDataPath(yamlPath, ev.benchmark_schema),

        // adata.obs keys
        sample_key:       keys.sample,
        cell_type_key:    keys.cell_type,
        batch_key:        keys.batch,
        batch_covariates: keys.batch_covariates,

        // metadata column lists
        samples_metadata_cols:       meta.samples_metadata_cols,
        accessible_metadata_columns: meta.accessible_metadata_columns,

        // preprocess settings
        sample_size_threshold: pp.sample_size_threshold,
        output_compression:    pp.output_compression,

        // evaluate settings
        root_sample:         ev.root_sample,
        trajectory_variable: ev.trajectory_variable,
        inverse_trajectory:  ev.inverse_trajectory,
        figure_format:       ev.figure_format,
    ]
}
