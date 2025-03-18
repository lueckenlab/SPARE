workflow run_wf {
  take:
    input_ch

  main:
    // ["hlca", [ input: "9f222629-9e39-47d0-b83f-e08d610c7479", source: "cxg", max_chunks: params.max_chunks ] ],
    // ["stephenson", [ input: "c7775e88-49bf-4ba2-a03b-93f00447c958", source: "cxg", max_chunks: params.max_chunks ] ],
    // ["onek1k", [ input: "3faad104-2ab8-4434-816d-474d8d2641db", source: "cxg", max_chunks: params.max_chunks ] ],

    combat_ch = 

      // Create a channel with two events
      // Each event contains a string (an identifier) and a file (input)
      Channel.fromList([
          ["combat", [ input: "https://zenodo.org/record/6120249/files/COMBAT-CITESeq-DATA.h5ad", source: "web", max_chunks: params.max_chunks ] ],
        ])

        // View channel contents
        | view { tup -> "Input: $tup" }
        | download_from_web.run(
            fromState: [ input: "input" ],
            toState: [ raw_dataset: "output" ]
        )
        | clean_combat.run(
          fromState: [ input: "raw_dataset" ],
          toState: [ cleaned_dataset: "output" ]
        )
        // TODO: add every dataset, process, represent, evaluate, then merge in one benchmarking table and return it here
          toState: { id, output, state ->
            output
          }
        )
        // TODO: add every dataset, process, represent, evaluate, then merge in one benchmarking table and return it here
        
        // download_from_cxg.run(
        //   runIf: { id, state ->
        //     state.source == "cxg"
        //   },
        //   fromState: { id, state ->
        //     def stateMapping = [
        //       "input": state.input,
        //       "output": "${id}/${id}.h5ad",
        //       "max_chunks": state.max_chunks,
        //     ]
        //     return stateMapping
        //   },
        //   toState: { id, output, state ->
        //     output
        //   }
        // )
        // // View channel contents
        // | view { tup -> "Output: $tup" }

  emit:
    combat_ch
      | map{ id, state -> [ "run", state ] }
}
