#!/usr/bin/env nextflow

/*
  ESVotron: A pipeline to quickly determine the true exact sequence variants from 16s rRNA amplicon reads.
*/

// Using DSL-2
nextflow.preview.dsl=2

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.help = false
params.output = './results'
params.manifest = null


// Function which prints help message text
def helpMessage() {
    log.info"""
    A workflow to obtain the exact sequence variants (ESV) of a community
    
    Usage:

    nextflow run composition_only.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results)
      --output_prefix       Text used as a prefix for summary HDF5 output files (default: geneshot)
      --nopreprocess        If specified, omit the preprocessing steps (removing adapters and human sequences)
      --savereads           If specified, save the preprocessed reads to the output folder (inside qc/)
      -w                    Working directory. Defaults to `./work`

    For preprocessing:
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
    
    Manifest file:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data for preprocessing is only accepted as paired reads.
      Reads are specified by columns, `read__1` and `read__2`.
      If index reads are provided, the column titles should be `I1` and `I2`
      If you wish to provide already processed data in fasta format, please include it in `R1` alone,
        with only *one* file specified per specimen.
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    output_folder = params.output.concat("/")
} else {
    output_folder = params.output
}

// Import the read_manifest module
include read_manifest from './modules/general'


workflow {
    main:

    // Phase 0: Validation of input data


}
