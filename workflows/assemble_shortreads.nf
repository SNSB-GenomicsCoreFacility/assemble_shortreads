/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                               } from '../modules/nf-core/fastqc/main'
include { FASTQC as ADAPTERFILTERFASTQC        } from '../modules/nf-core/fastqc/main'
include { FASTQC as DEDUPFASTQC                } from '../modules/nf-core/fastqc/main'
include { FASTQC as REPAIREDFASTQC             } from '../modules/nf-core/fastqc/main'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { BBMAP_BBDUK                          } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_REPAIR                         } from '../modules/nf-core/bbmap/repair/main'
include { FASTP                                } from '../modules/nf-core/fastp/main'
include { PARDRE                               } from '../modules/local/pardre/main'
include { SPADES                               } from '../modules/nf-core/spades/main'
include { QUAST                                } from '../modules/nf-core/quast/main'
include { MINIMAP2_ALIGN                       } from '../modules/nf-core/minimap2/align/main'
include { MUMMER                               } from '../modules/nf-core/mummer/main'
include { BWAMEM2_INDEX                        } from '../modules/nf-core/bwamem2/index/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_assemble_shortreads_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLE_SHORTREADS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // read adapter fasta files 

    fasta_adapter=Channel.fromPath(params.reference_adapter_fas,checkIfExists: true)
    
    // prepare input channel for bbmap adapter trimming

    prech_adapter_files = ch_samplesheet.combine(fasta_adapter)


    prech_adapter_files.multiMap{meta,reads,fasta ->
        //first_ch: tuple(meta, reads)
        first_ch: tuple([id:meta.id+"_adapter_filtered"]+[single_end:false], reads)
        second_ch: fasta
        }.set{ ch_adapter_files }

    //ch_adapter_files.first_ch.view()
    if(params.filtering_tool == "bbduk"){
        //
        // MODULE: BBMAP_BBDUK
        //
        BBMAP_BBDUK(
            ch_adapter_files.first_ch,
            ch_adapter_files.second_ch
         )
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{it[1]})
        ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        ch_adapter_fastqc = BBMAP_BBDUK.out.reads
    }

    else{
        //
        // MODULE: FASTP
        //
        FASTP(
            ch_adapter_files.first_ch,
            ch_adapter_files.second_ch,
            Channel.value(false),
            Channel.value(false),
            Channel.value(false)
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_adapter_fastqc = FASTP.out.reads
        }

        //
        // MODULE: ADAPTER_FILTER_FASTQC
        //
        ADAPTERFILTERFASTQC(
            ch_adapter_fastqc
        )

        ch_multiqc_files = ch_multiqc_files.mix(ADAPTERFILTERFASTQC.out.zip.collect{it[1]})

        //
        // MODULE : PARDRE
        //
        PARDRE(
            ch_adapter_fastqc
        )
        
        ch_versions = ch_versions.mix(PARDRE.out.versions.first())

        prech_pardre_reads = PARDRE.out.reads

        prech_pardre_reads.map{
            meta,reads ->tuple([id:meta.id+"_dedup_pardre"]+[single_end:false], reads)
        }.set{ ch_pardre_reads }

        //
        //MODULE: DEDUPFASTQC
        //
        DEDUPFASTQC(
            ch_pardre_reads
        )

        ch_multiqc_files = ch_multiqc_files.mix(DEDUPFASTQC.out.zip.collect{it[1]})
    

        //
        //MODULE: BBMAP_REPAIR
        //
        BBMAP_REPAIR(
            ch_pardre_reads,
            channel.value(false)
        )

        BBMAP_REPAIR.out.repaired.map{
            meta,reads->tuple([id:meta.id+"_repaired"]+[single_end:false],reads)
        }.set{ch_bbmap_repaired_reads}

        //
        //MODULE: REPAIREDFASTQC
        //
        REPAIREDFASTQC(
            ch_bbmap_repaired_reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(REPAIREDFASTQC.out.zip.collect{it[1]})

        //
        //MODULE: SPADES
        //
        SPADES(
            ch_bbmap_repaired_reads.map{meta,reads->tuple(meta,reads,[],[])},
            [],
            []
        )

        if(params.reference_fasta){
            if(params.reference_fasta.endsWith(".map")){
                    Channel
                        .fromPath(params.reference_fasta)
                        .splitCsv(header:true,sep:",")
                        .map{row ->
                            tuple([id:row.sample],row.reference)
                        }
                        .set{meta_reference}

                    SPADES.out.scaffolds
                        .map{meta,scaffolds ->
                            def new_id = meta.id.replaceFirst(/_adapter_filtered_dedup_pardre_repaired$/,'')
                            def new_meta = [:]
                            new_meta.id = new_id
                            [new_meta,scaffolds]
                        }
                        .set{nmeta_scaffolds}
                    
                    nmeta_scaffolds
                        .join(meta_reference)
                        .multiMap{meta,scaffolds,reference ->
                                 first_ch: tuple(meta,scaffolds)
                                 second_ch: tuple(meta,reference)
                        }
                        .set{ch_minimap2_align}

            }

            else {
            reference_fasta = Channel.fromPath(params.reference_fasta)
            prech_reference_fasta = reference_fasta.map{fasta->tuple([id:"reference"],fasta)}

            prech_scaffolds_fasta = SPADES.out.scaffolds.map{meta, scaffolds -> tuple([id:meta.id]+[single_end:"false"],scaffolds)}

            //MODULE: MINIMAP2_ALIGN


                prech_minimap2_align = prech_scaffolds_fasta.combine(prech_reference_fasta)

                prech_minimap2_align.multiMap{ meta, scaffolds, meta2, ref -> 
                    first_ch: tuple(meta,scaffolds)
                    second_ch: tuple(meta,ref)
                }.set{ch_minimap2_align}
            }

            if(params.galignment_tool == "minimap2" ){

                //
                //MODULE: MINIMAP2_ALIGN
                //

                MINIMAP2_ALIGN(
                    ch_minimap2_align.first_ch,
                    ch_minimap2_align.second_ch,
                    channel.value(false),
                    channel.value("bai"),
                    channel.value(true),
                    channel.value(true)
                )
           }

           if(params.galignment_tool == "mummer" ){
                   
                   //
                   // MODULE: MUMMER
                   //

                ch_minimap2_align.first_ch
                    .join(ch_minimap2_align.second_ch)
                    .set{ ch_mummer }
                

                MUMMER(
                    ch_mummer
                )

            if(params.consensus_fasta==true){

                if(params.reference_fasta.endsWith(".map")){
                    // Read the CSV and parse rows into maps
                    Channel
                        .fromPath(params.reference_fasta)
                        .splitCsv(header:true)
                        .map { row -> 
                            [ 
                                sample: row.sample, 
                                reference: row.reference, 
                                reference_idx: row.reference_idx 
                            ]
                        }
                        .branch{it -> 
                            withIdx: it.reference_idx != "none"  
                            withoutIdx: it.reference_idx == "none"
                        }.set{ref}

                    cp_ref_without_idx = ref.withoutIdx.map{it->tuple([id:it.sample],it.reference)}
                    cp_ref_with_idx = ref.withIdx.map{it->tuple([id:it.sample],it.reference,it.reference_idx)}


                    //
                    //BWAMEM2_INDEX
                    //

                    BWAMEM2_INDEX(
                        cp_ref_without_idx
                    )

                    prech_ref_with_idx = cp_ref_without_idx.join(BWAMEM2_INDEX.out.index)

                }
            }

            }

        
        }

        prech_quast = SPADES.out.scaffolds.map{meta,assembly->assembly}.collect()

        ch_quast = prech_quast.map{assemblies->tuple([id:"consensus"],assemblies)}

        QUAST(
            ch_quast,
            [[],[]],
            [[],[]]
        )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'assemble_shortreads_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
