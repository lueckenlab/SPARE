include { expandDatasetInfo } from "./dataset_info.nf"

// Fill any missing arg from state.dataset_info; explicit state values win.
def applyDatasetInfo = { id, state ->
    if (!state.dataset_info) return [id, state]
    def info = readYaml(state.dataset_info)
    def expanded = expandDatasetInfo(info, state.dataset_info as String)
    def extras = [
        input:                       expanded.representations_path,
        metadata_path:               expanded.metadata_path,
        cell_type_key:               expanded.cell_type_key,
        accessible_metadata_columns: expanded.accessible_metadata_columns,
        benchmark_schema:            expanded.benchmark_schema,
        root_sample:                 expanded.root_sample,
        trajectory_variable:         expanded.trajectory_variable,
        inverse_trajectory:          expanded.inverse_trajectory,
        figure_format:               expanded.figure_format,
    ].findAll { it.value != null }
    return [id, extras + state.findAll { it.value != null }]
}

workflow run_wf {
    take:
        input_ch

    main:
        expanded_ch = input_ch
            | map(applyDatasetInfo)
            | view { id, state ->
                "\n[evaluate] ID: $id\n  input=${state.input}\n  method_outputs=${state.method_outputs}\n  benchmark_schema=${state.benchmark_schema}\n  trajectory_variable=${state.trajectory_variable}"
            }

        // Aggregate first if per-method CSVs were supplied.
        aggregated_ch = expanded_ch
            | branch { id, state ->
                aggregate: state.method_outputs != null && state.method_outputs.size() > 0
                skip:      true
            }

        aggregated_ch.aggregate
            | aggregate_representations.run(
                fromState: { id, state ->
                    [
                        input: state.method_outputs,
                        output: state.input,
                        metadata_path: state.metadata_path,
                        cell_type_key: state.cell_type_key,
                        accessible_metadata_columns: state.accessible_metadata_columns,
                    ]
                },
                toState: { id, output, state ->
                    state + [input: output.output]
                }
            )
            | set { post_aggregate_ch }

        ready_ch = post_aggregate_ch.mix(aggregated_ch.skip)

        output_ch = ready_ch
            | evaluate.run(
                fromState: { id, state ->
                    def args = [
                        input: state.input,
                        output_dir: state.output_dir,
                        benchmark_schema: state.benchmark_schema,
                        figure_format: state.figure_format ?: "png",
                        inverse_trajectory: state.inverse_trajectory ?: false,
                    ]
                    if (state.root_sample != null)         args.root_sample = state.root_sample
                    if (state.trajectory_variable != null) args.trajectory_variable = state.trajectory_variable
                    return args
                },
                toState: { id, output, state ->
                    state + [output_dir: output.output_dir]
                }
            )

    emit:
        output_ch
}
