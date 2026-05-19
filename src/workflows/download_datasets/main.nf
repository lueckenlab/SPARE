include { expandDatasetInfo } from "./dataset_info.nf"

// If state.dataset_info is set, fill missing fields from the YAML.
// Individual args (input/source/output) win where present.
def applyDatasetInfo = { id, state ->
    if (!state.dataset_info) return [id, state]
    def info = readYaml(state.dataset_info)
    def expanded = expandDatasetInfo(info, state.dataset_info as String)
    def src = info.source ?: [:]
    def extras = [
        input:  src.url ?: src.cxg_id,
        source: src.type,
        output: expanded.raw_path,
    ].findAll { it.value != null }
    return [id, extras + state.findAll { it.value != null }]
}

workflow run_wf {
  take:
    input_ch

  main:

    output_ch =
      input_ch
        | map(applyDatasetInfo)
        | view { tup -> "Input: $tup" }
        | download_from_web.run(
          runIf: { id, state ->
            state.source == "web" && !file(state.output).exists()
          },
          fromState: { id, state ->
            [input: state.input, output: state.output, max_chunks: state.max_chunks]
          },
          toState: { id, output, state ->
            output
          }
        )
        | download_from_cxg.run(
          runIf: { id, state ->
            state.source == "cxg" && !file(state.output).exists()
          },
          fromState: { id, state ->
            [input: state.input, output: state.output, max_chunks: state.max_chunks]
          },
          toState: { id, output, state ->
            output
          }
        )
        | view { tup -> "Output: $tup" }

  emit:
    output_ch
}
