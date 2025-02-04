workflow run_wf {
  take:
    input_ch

  main:

    output_ch = 
      input_ch
        // View channel contents
        | view { tup -> "Input: $tup" }
        | download_from_web.run(
          runIf: { id, state ->
            state.source == "web"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "max_chunks": state.max_chunks,
            ]
            return stateMapping
          },
          toState: { id, output, state ->
            output
          }
        )
        | download_from_cxg.run(
          runIf: { id, state ->
            state.source == "cxg"
          },
          fromState: { id, state ->
            def stateMapping = [
              "input": state.input,
              "output": state.output,
              "max_chunks": state.max_chunks,
            ]
            return stateMapping
          },
          toState: { id, output, state ->
            output
          }
        )
        // View channel contents
        | view { tup -> "Output: $tup" }

  emit:
    output_ch
}
