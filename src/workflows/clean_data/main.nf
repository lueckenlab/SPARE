workflow run_wf {
  take:
    input_ch

  main:

    output_ch = 
      input_ch
        // View channel contents
        | view { tup -> "Input: $tup" }
        | clean_combat.run(
          runIf: { id, state ->
            id == "combat"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "output_compression": state.output_compression,
            ]
            return stateMapping
          },
          toState: { id, output, state ->
            output
          }
        )
        | clean_stephenson.run(
          runIf: { id, state ->
            id == "stephenson"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "output_compression": state.output_compression,
            ]
            return stateMapping
          },
        )
        | clean_hlca.run(
          runIf: { id, state ->
            id == "hlca"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "output_compression": state.output_compression,
            ]
            return stateMapping
          },
        )
        | clean_onek1k.run(
          runIf: { id, state ->
            id == "onek1k"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "output_compression": state.output_compression,
            ]
            return stateMapping
          },
        )
        // View channel contents
        | view { tup -> "Output: $tup" }

  emit:
    output_ch
}
