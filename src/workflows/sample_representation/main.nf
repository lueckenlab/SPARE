def methodMap = [
    pseudobulk: pseudobulk,
    random_vector: random_vector
]

def runTriplets() {
    workflow runTripletsWf {
        take: input_ch
        main:
            method_ch = input_ch.flatMap { id, state ->
                state.triplets.collect { triplet ->
                    def method_name = triplet.method.config.name
                    def experiment_name = triplet.id.tokenize('.')[-1]
                    [ 
                        "${id}.${method_name}.${experiment_name}",
                        triplet.state,
                        [ 
                            output: triplet.state.output,
                            aggregated_output: triplet.state.aggregated_output,
                            metadata_path: triplet.state.metadata_path,
                            _meta: [join_id: id]
                        ],
                        triplet.method
                    ]
                }
            }
            
            // Use flatMap to handle method execution channels
            output_ch = method_ch.flatMap { id, state, orig_state, method ->
                println "Running method: ${method.config?.name} with ID: ${id}"
                
                // Execute the method and get its OUTPUT CHANNEL
                def result_ch = method.run(
                    auto: [simplifyInput: false, simplifyOutput: false],
                    args: state
                )
                
                // Transform the method's output into the required format
                result_ch.map { method_output ->
                    [ 
                        orig_state._meta.join_id, 
                        [ 
                            output: method_output.output, // Use actual output from method
                            aggregated_output: orig_state.aggregated_output,
                            metadata_path: orig_state.metadata_path,
                            _meta: orig_state._meta
                        ] 
                    ]
                }
            }
            | groupTuple()
            | map { id, states ->
                [ 
                    id, 
                    [ 
                        output: states.collect { it.output },
                        aggregated_output: states[0].aggregated_output,
                        metadata_path: states[0].metadata_path,
                        _meta: [join_id: id]
                    ] 
                ]
            }
            
        emit:
            output_ch
    }
    return runTripletsWf
}

workflow run_wf {
    take:
        input_ch

    main:
        output_ch = input_ch
            | view { id, state -> 
                "\n[Initial Input] ID: $id\nState: $state" 
            }
            
            // Read method parameters and convert to triplets
            | map { id, state ->
                def method_params_file = file(state.method_params)
                def method_params = method_params_file.exists() ? readYaml(method_params_file) : [:]
                
                // Convert nested map to list of triplets [method_name, experiment_name, params]
                def method_triplets = []
                method_params.each { method_name, experiments ->
                    // Get the method component from our map
                    def method = methodMap[method_name]
                    if (!method) {
                        error "Method $method_name not found in available methods"
                    }
                    
                    experiments.each { experiment_name, params ->
                        method_triplets << [
                            method: method,
                            id: id,
                            state: [input: state.input] + params + [
                                output: state.output.replace('.h5ad', "_${method_name}_${experiment_name}.csv"),
                                aggregated_output: state.output,
                                sample_key: state.sample_key,
                                cell_type_key: state.cell_type_key,
                                metadata_path: state.metadata_path
                            ]
                        ]
                    }
                }
                
                [id, state + [triplets: method_triplets]]
            }
            | view { id, state -> 
                "\n[After Reading Files] ID: $id\nTriplets: ${state.triplets}" 
            }
            // Add metadata
            | map { id, state ->
                [id, state + ["_meta": [join_id: id]]]
            }
            | view { id, state -> 
                "\n[After Adding Metadata] ID: $id\nMeta: ${state._meta}" 
            }
            // Run all triplets
            | runTriplets()
            | view { id, state -> 
                "\n[Final Output] ID: $id\nOutput: ${state.output}, state: ${state}" 
            }
            | aggregate_representations.run(
                fromState: { id, state ->
                    [
                        input: state.output,  // Individual method output CSV
                        output: state.aggregated_output,  // Final aggregated h5ad
                        metadata_path: state.metadata_path,
                        cell_type_key: state.cell_type_key
                    ]
                }
            )

    emit:
        output_ch
}
