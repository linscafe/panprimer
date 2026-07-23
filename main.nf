#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Pangenome PCR primer design — one process per pipeline stage (see the plan).
// Per-haplotype stages (PROJECT, EVALUATE) fan out over the haplotype channel; MASK_DESIGN
// and AGGREGATE are the fan-in points. Each process shells out to a `pangenome-primer`
// subcommand with JSON artifacts as the inter-stage contract.

process ANCHOR {
    input:
      val target
      path chm13
      path cfg
    output:
      path 'anchor.json'
    script:
    """
    pangenome-primer anchor --target '${target}' --chm13 ${chm13} --config ${cfg} --out anchor.json
    """
}

process PROJECT {
    tag "${hap_id}"
    input:
      path anchor
      tuple val(hap_id), path(hap_fasta)
    output:
      tuple val(hap_id), path("proj.${hap_id.replace('#','_')}.json"), path(hap_fasta)
    script:
    """
    pangenome-primer project --anchor ${anchor} --hap-fasta ${hap_fasta} \
        --hap-id '${hap_id}' --out proj.${hap_id.replace('#','_')}.json
    """
}

process MASK_DESIGN {
    input:
      path anchor
      path projections
    output:
      path 'candidates.json'
    script:
    def proj_args = projections.collect { "--projection ${it}" }.join(' ')
    """
    pangenome-primer design --anchor ${anchor} ${proj_args} --out candidates.json
    """
}

process EVALUATE {
    tag "${hap_id}"
    input:
      path candidates
      tuple val(hap_id), path(proj), path(hap_fasta)
      val mode
    output:
      path "result.${hap_id.replace('#','_')}.json"
    script:
    def mode_arg = mode ? "--mode ${mode}" : ''
    """
    pangenome-primer evaluate --candidates ${candidates} --projection ${proj} \
        --hap-fasta ${hap_fasta} ${mode_arg} --out result.${hap_id.replace('#','_')}.json
    """
}

process AGGREGATE {
    publishDir "${params.outdir}", mode: 'copy'
    input:
      path anchor
      path candidates
      path results
    output:
      path 'report/*'
    script:
    def res_args = results.collect { "--result ${it}" }.join(' ')
    """
    pangenome-primer aggregate --anchor ${anchor} --candidates ${candidates} \
        ${res_args} --outdir report
    """
}

workflow {
    if( !params.target || !params.chm13 )
        error "Required: --target 'chr:start-end'|FASTA and --chm13 <CHM13 FASTA>"

    cfg   = file(params.config)
    chm13 = file(params.chm13)

    // haplotype channel: (hap_id, fasta) from samples.tsv (url column = local path)
    haplos = Channel.fromPath(params.samples)
        .splitCsv(sep: '\t', header: false)
        .filter { it && !it[0].startsWith('#') && it[0] != 'sample' }
        .map { row -> tuple("${row[0]}#hap${row[1]}", file(row[3])) }

    anchor      = ANCHOR(params.target, chm13, cfg)
    projections = PROJECT(anchor, haplos)                 // (hap_id, proj.json, hap_fasta)
    proj_jsons  = projections.map { it[1] }.collect()
    candidates  = MASK_DESIGN(anchor, proj_jsons)
    results     = EVALUATE(candidates, projections, params.mode ?: '')
    AGGREGATE(anchor, candidates, results.collect())
}
