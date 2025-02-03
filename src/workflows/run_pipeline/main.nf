workflow run_wf {
  take:
    input_ch

  main:
    input_ch | view { "Input channel contains: $it" }

    // Create a dummy input file. Otherwise the pipeline doesn't run even though the input is not required for "download" components
    dummy_file = file("${projectDir}/NO_FILE")
    dummy_file.text = ""  // Create an empty file

    combat_ch = Channel.fromList([
      ["run", [input: dummy_file]]  // Changed ID to "run" to match input channel
    ]) | download_combat.run(
      fromState: [ output_dir: "data" ]
    )

    stephenson_ch = Channel.fromList([
      ["run", [input: dummy_file]]  // Changed ID to "run" to match input channel
    ]) | download_stephenson.run(
      fromState: [ output_dir: "data" ]
    )

    // Combine both outputs
    output_ch = combat_ch.mix(stephenson_ch)

    output_ch | view { "Output channel contains: $it" }

  emit:
    output_ch
}
