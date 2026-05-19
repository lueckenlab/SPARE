include { expandDatasetInfo } from "./dataset_info.nf"

def applyDatasetInfo = { id, state ->
    if (!state.dataset_info) return [id, state]
    def info = readYaml(state.dataset_info)
    def expanded = expandDatasetInfo(info, state.dataset_info as String)
    def extras = [
        input:  expanded.raw_path,
        output: expanded.cleaned_path,
        output_compression: expanded.output_compression,
    ].findAll { it.value != null }
    return [id, extras + state.findAll { it.value != null }]
}

def cleanArgs = { id, state ->
    [
        input:  state.input,
        output: state.output,
        output_compression: state.output_compression,
    ]
}

workflow run_wf {
  take:
    input_ch

  main:

    output_ch =
      input_ch
        | map(applyDatasetInfo)
        | view { tup -> "Input: $tup" }
        | clean_combat.run(
          runIf:    { id, state -> id == "combat" },
          fromState: cleanArgs,
        )
        | clean_stephenson.run(
          runIf:    { id, state -> id == "stephenson" },
          fromState: cleanArgs,
        )
        | clean_hlca.run(
          runIf:    { id, state -> id == "hlca" },
          fromState: cleanArgs,
        )
        | clean_onek1k.run(
          runIf:    { id, state -> id == "onek1k" },
          fromState: cleanArgs,
        )
        | clean_ticatlas.run(
          runIf:    { id, state -> id == "ticatlas" },
          fromState: cleanArgs,
        )
        | clean_lupus.run(
          runIf:    { id, state -> id == "lupus" },
          fromState: cleanArgs,
        )
        | clean_imm_of_aging.run(
          runIf:    { id, state -> id == "imm_of_aging" },
          fromState: cleanArgs,
        )
        | clean_sound_life.run(
          runIf:    { id, state -> id == "sound_life" },
          fromState: cleanArgs,
        )
        | view { tup -> "Output: $tup" }

  emit:
    output_ch
}
