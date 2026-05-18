include { expandDatasetInfo } from "${params.rootDir}/src/workflows/utils/dataset_info.nf"

methods = [
    pseudobulk,
    random_vector,
    mrvi,
    scpoli,
    grouped_pseudobulk,
    pilot,
    mofa,
    gloscope,
    cell_group_composition,
]

// If state.dataset_info is set, fill missing fields from the YAML.
// Existing state values win when they are non-null.
def applyDatasetInfo = { id, state ->
    if (!state.dataset_info) return [id, state]
    def info = readYaml(state.dataset_info)
    def expanded = expandDatasetInfo(info, state.dataset_info as String)
    def extras = [
        input:                       expanded.processed_path,
        output:                      expanded.representations_path,
        metadata_path:               expanded.metadata_path,
        sample_key:                  expanded.sample_key,
        cell_type_key:               expanded.cell_type_key,
        accessible_metadata_columns: expanded.accessible_metadata_columns,
    ].findAll { it.value != null }
    return [id, extras + state.findAll { it.value != null }]
}

workflow run_wf {
    take:
        input_ch

    main:
        // First expand the experiments
        expanded_ch = input_ch
            | map(applyDatasetInfo)
            | view { id, state ->
                "\n[Initial Input] ID: $id\nState: $state"
            }
            | map { id, state ->
                def experiments = readYaml(state.method_params)
                [id, state + ["_meta": [join_id: id], "experiments": experiments]]
            }
            | view { id, state ->
                "\n[After Adding Metadata] ID: $id\nMeta: ${state._meta}"
            }
            | flatMap { id, state ->
                def runs = []
                methods.each { method ->
                    def methodExperiments = state.experiments[method.config.name]
                    methodExperiments.each { experimentName, params ->
                        def new_state = state + [
                            "current_method": method.config.name,
                            "current_experiment": experimentName,
                            "current_params": params,
                            "_meta": state._meta + [original_id: id]
                        ]
                        def run_id = "${id}.${method.config.name}.${experimentName}".toString()
                        runs << [
                            run_id,
                            new_state
                        ]
                    }
                }
                return runs
            }
            | view { id, state ->
                "\n[Expanded Experiments] ID: $id\nMethod: ${state.current_method}\nExperiment: ${state.current_experiment}"
            }

        // Run methods on expanded experiments
        method_outputs_ch = expanded_ch
            | runEach(
                components: methods,
                filter: { id, state, comp ->
                    comp.config.name == state.current_method
                },
                id: { id, state, comp ->
                    id.toString()
                },
                fromState: { id, state, comp ->
                    def methodName = state.current_method
                    def experimentName = state.current_experiment
                    def experimentParams = state.current_params

                    def new_args = [
                        input: state.input,
                        metadata_path: state.metadata_path,
                        cell_type_key: state.cell_type_key,
                        sample_key: state.sample_key,
                        output: "${methodName}_${experimentName}.csv",
                    ]
                    new_args.putAll(experimentParams)
                    return new_args
                },
                toState: { id, output, state, comp ->
                    state + [
                        method_id: "${state.current_method}_${state.current_experiment}".toString(),
                        method_output: output.output
                    ]
                }
            )
            | view { id, state ->
                "\n[Method Output] ID: $id\nOutput: ${state.method_output}"
            }

        // Group and aggregate results
        output_ch = method_outputs_ch
            | map { id, state ->
                [state._meta.original_id.toString(), [id, state]]
            }
            | groupTuple()
            | map { original_id, grouped_outputs ->
                def first_state = grouped_outputs[0][1]
                [
                    original_id,
                    [
                        method_outputs: grouped_outputs.collect { it[1].method_output },
                        output: first_state.output,
                        metadata_path: first_state.metadata_path,
                        cell_type_key: first_state.cell_type_key,
                        accessible_metadata_columns: first_state.accessible_metadata_columns,
                        _meta: [join_id: original_id]
                    ]
                ]
            }
            | view { id, state ->
                "\n[Grouped Outputs] ID: $id\nOutputs: ${state.method_outputs}"
            }
            | aggregate_representations.run(
                fromState: { id, state ->
                    [
                        input: state.method_outputs,
                        output: state.output,
                        metadata_path: state.metadata_path,
                        cell_type_key: state.cell_type_key,
                        accessible_metadata_columns: state.accessible_metadata_columns
                    ]
                }
            )

    emit:
        output_ch
}
