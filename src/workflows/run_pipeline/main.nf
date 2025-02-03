workflow run_wf {
  take:
    input_ch

  main:
    input_ch | view { "Input channel contains: $it" }

    // Create a dummy input file. Otherwise the pipeline doesn't run even though the input is not required for "download" components
    dummy_file = file("${projectDir}/NO_FILE")
    dummy_file.text = ""  // Create an empty file

    combat_ch = Channel.fromList([
      ["run", [input: dummy_file]]
    ]) | download_combat.run(
      fromState: [ output: "${params.output_dir}/combat/combat.h5ad" ]
    ) | view { "COMBAT channel contains: $it" }

    stephenson_ch = Channel.fromList([
      ["run", [input: dummy_file]]
    ]) | download_stephenson.run(
      fromState: [ output: "${params.output_dir}/stephenson/stephenson.h5ad" ]
    ) | view { "Stephenson channel contains: $it" }

  emit:
    combat_ch
    stephenson_ch
}
