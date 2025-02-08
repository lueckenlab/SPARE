methods = [
    pseudobulk, 
    random_vector,
    mrvi,
    scpoli,
    grouped_pseudobulk,
]

workflow run_wf {
    take:
        input_ch

    main:
        // First expand the experiments
        expanded_ch = input_ch
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
                        // Store original ID in metadata for later grouping
                        def new_state = state + [
                            "current_method": method.config.name,
                            "current_experiment": experimentName,
                            "current_params": params,
                            "_meta": state._meta + [original_id: id]
                        ]
                        // Convert GString to plain String using toString()
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
                // Only run a method if it matches the current experiment's method
                filter: { id, state, comp ->
                    comp.config.name == state.current_method
                },
                id: { id, state, comp ->
                    id.toString() // Ensure we return a plain String
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
                    
                    // Add experiment-specific parameters
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
                // Extract original ID for grouping
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
                        obs_columns: first_state.obs_columns,
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
                        obs_columns: state.obs_columns
                    ]
                }
            )

    emit:
        output_ch
}
