#!/usr/bin/env nextflow

/*
  ESVotron: A pipeline to quickly determine the true exact sequence variants from 16s rRNA amplicon reads.
*/

// Using DSL-2
nextflow.preview.dsl=2

// Containers
container__vsearch = "golob/vsearch:2.7.1_bcw_0.2.0"
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"
container__fastcombineseqtab = "golob/dada2-fast-combineseqtab:0.5.0__1.12.0__BCW_0.3.1"

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.help = false
params.output = './results'
params.manifest = null

// Operational params
params.collector_iteration_cutoff = 0.0001
params.collector_min_reads = 500

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
      -w                    Working directory. Defaults to `./work`
    
    Manifest file:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data for preprocessing is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
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
include countReads from './modules/general'


process pairFastq {
    container "${container__vsearch}"
    label = 'multithread'
    publishDir "${params.output}/quickmerged_reads/", mode: 'copy'

    input:
        tuple specimen, file(R1), file(R2)
    
    output:
        tuple specimen, file("${specimen}.fasta.gz")

    """
    set -e 

    vsearch \
    --threads=${task.cpus} \
    --fastq_mergepairs \
    ${R1} --reverse ${R2} \
    --fastaout ${specimen}.fasta \
    --fastq_maxns 0 \
    --fastq_maxdiffs 0

    gzip ${specimen}.fasta
    """
}

process find_saturation_read {
    container "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/collector/", mode: 'copy'

    input:
        tuple specimen, file(Reads)
    
    output:
        tuple specimen, stdout, file(Reads), file("${specimen}.esv_counts.csv"), file("${specimen}.collector.csv")

"""
#!/usr/bin/env python3

from fastalite import fastalite
import gzip
import csv
from collections import defaultdict

with gzip.open('${Reads}', 'rt') as fasta_h:
    esv_count = defaultdict(int)
    rareifaction_curve = []
    saturation_reads = None
    for sr_i, sr in enumerate(fastalite(fasta_h)):
        esv_count[sr.seq] += 1
        total_reads = sr_i + 1
        esv_total = len(esv_count)
        esv_singleton = len([
            v for v in 
            esv_count.values()
            if v == 1
        ])
        rareifaction_curve.append([
            total_reads,
            esv_total,
            esv_singleton,
            1 - esv_singleton / total_reads
        ])
        if saturation_reads is None and total_reads > ${params.collector_min_reads} and abs(rareifaction_curve[-1][3] - rareifaction_curve[-2][3]) <= ${params.collector_iteration_cutoff}:
            saturation_reads = total_reads
with open('${specimen}.collector.csv', 'wt') as collect_h:
    collect_w =  csv.writer(collect_h)
    collect_w.writerow([
        'reads_total',
        'esv_total',
        'esv_singleton',
        'goods',
        'threshold'
    ])
    for r in rareifaction_curve:
        collect_w.writerow([
            r[0],
            r[1],
            r[2],
            r[3],
            saturation_reads
        ])

with open('${specimen}.esv_counts.csv', 'wt') as esv_count_h:
    esv_w =  csv.writer(esv_count_h)
    esv_w.writerow([
        '',
        '${specimen}',
    ])
    for esv, count in esv_count.items():
        esv_w.writerow([
            esv,
            count,
        ])

print(saturation_reads)
"""
}

process filter_counts_by_threshold {
    container "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/filtered_counts/", mode: 'copy'

    input:
        tuple specimen, val(threshold), file(Reads), file(esv_count), file(collector)
    
    output:
        tuple specimen, file("${specimen}.esv_count_filtered.csv")
        

"""
#!/usr/bin/env python3

import csv
with open('${specimen}.esv_count_filtered.csv', 'wt') as out_h:
    out_w = csv.writer(out_h)
    esv_count_r = csv.reader(
        open('${esv_count}', 'rt')
    )
    header = next(esv_count_r)
    out_w.writerow(header)
    out_w.writerows(
        [
            r for r in esv_count_r
            if int(r[1]) >= (${threshold} / 3)
        ]
    )


"""
}

process combine_esv_counts {
    container "${container__fastcombineseqtab}"
    label = 'io_limited'
    publishDir "${params.output}/", mode: 'copy'

    input:
        file(esv_counts)
    
    output:
        file("seqtab.csv")
        

"""
#!/usr/bin/env python3

import pandas as pd
import csv
import random

esv_count_files = "${esv_counts}".split()
random.shuffle(esv_count_files)
total_counts_d = {}
for ecf in esv_count_files:
    ecf_r = csv.reader(open(ecf, 'rt'))
    specimen = next(ecf_r)[1]
    total_counts_d[specimen] = {
        r[0]: r[1]
        for r in ecf_r
    }

seqtab = pd.DataFrame(total_counts_d).fillna(0).astype(int)

# Get rid of empty specimens

seqtab = seqtab[
    seqtab.columns[seqtab.sum() > 0]
]
# Get rid of global singletons
seqtab = seqtab.loc[(seqtab.T > 0).sum() > 0]
seqtab.to_csv('seqtab.csv')

"""
}

workflow {
    main:

    // Phase 0: Validation of input data
    manifest = read_manifest(
        Channel.from(
            file(params.manifest)
        )
    )

    manifest.valid_paired_indexed.mix(manifest.valid_paired)
        .map{ r -> [r.specimen, file(r.R1), file(r.R2)] }
        .set { manifest_valid_ch }

    // Phase 1: Quick n dirty collector's curves
    specimen_merged_reads = pairFastq(manifest_valid_ch)

    saturation_thresholds = find_saturation_read(specimen_merged_reads)
        .map{ 
            r -> [
                r[0],
                r[1].strip(),
                r[2],
                r[3],
                r[4]
            ]
        }
        .branch {
            failed:  it[1] == 'None'
            converged: true
        }
    
    saturation_thresholds_passed = saturation_thresholds.converged
        .map{ r-> [
            r[0],
            r[1].toInteger(),
            file(r[2]),
            file(r[3]),
            file(r[4])
        ]}
    
    filtered_counts = filter_counts_by_threshold(saturation_thresholds_passed)
    
    combine_esv_counts(
            filtered_counts
        .map { r -> file(r[1]) }
        .toList()
    )


    
}
