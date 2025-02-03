workflow run_wf {
  take:
    input_ch

  main:
    input_ch | view { "Input channel contains: $it" }

    // Create a dummy input file. Otherwise the pipeline doesn't run even though the input is not required for "download" components
    dummy_file = file("${projectDir}/NO_FILE")
    dummy_file.text = ""  // Create an empty file

    // Create input channels for each download component
    combat_input = Channel.fromList([
      ["run", [input: dummy_file]]
    ])

    stephenson_input = Channel.fromList([
      ["run", [input: dummy_file]]
    ])

    // Run the download components with proper input channels
    combat_ch = combat_input | 
      download_combat.run(
        fromState: [ output: "${params.output_dir}/combat/combat.h5ad" ]
      ) | 
      view { "COMBAT channel contains: $it" }

    stephenson_ch = stephenson_input | 
      download_stephenson.run(
        fromState: [ output: "${params.output_dir}/stephenson/stephenson.h5ad" ]
      ) | 
      view { "Stephenson channel contains: $it" }

  emit:
    combat_ch
    stephenson_ch
}
