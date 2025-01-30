nextflow.enable.dsl = 2

workflow download_all {
    take:
        input_ch

    main:
        combat_ch = input_ch \
            | download_combat.run(
                fromState: [ output: "${params.output_dir}/combat.h5ad", compression: params.compression ]
            )

        onek1k_ch = input_ch \
            | download_onek1k.run(
                fromState: [ output: "${params.output_dir}/onek1k.h5ad", compression: params.compression ]
            )

        stephenson_ch = input_ch \
            | download_stephenson.run(
                fromState: [ output: "${params.output_dir}/stephenson.h5ad", compression: params.compression ]
            )

        hlca_ch = input_ch \
            | download_hlca.run(
                fromState: [ output: "${params.output_dir}/hlca.h5ad", compression: params.compression ]
            )

        output_ch = combat_ch.mix(onek1k_ch, stephenson_ch, hlca_ch)

    emit:
        output_ch
} 
