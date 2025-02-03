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
      fromState: [ output: "data/combat/combat.h5ad" ]
    )

    stephenson_ch = Channel.fromList([
      ["run", [input: dummy_file]]
    ]) | download_stephenson.run(
      fromState: [ output: "data/stephenson/stephenson.h5ad" ]
    )

    // Combine both outputs
    output_ch = combat_ch.mix(stephenson_ch)

    output_ch | view { "Output channel contains: $it" }

  emit:
    output_ch
}
