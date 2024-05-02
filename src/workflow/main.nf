// process CreateDataDirectory {
//     output:
//     path "data"

//     script:
//     """
//     mkdir -p data
//     """
// }
params.dataset = "${params.rootDir}/data/combat_processed.h5ad"

workflow run_wf {
  take:
    input_ch

  main:
    // CreateDataDirectory()
    output_ch =

      input_ch
 
        | preprocess_combat.run(
            fromState: [ dataset: "dataset" ],
            // directives: [ publishDir: "output/"]
            // toState: { id, result, state -> state + result }
          )

        | represent_combat.run(
            fromState: [ dataset: "dataset" ],
            // directives: [ publishDir: "output/"]  
            // toState: { id, result, state -> result }
          )

  emit:
    output_ch
}