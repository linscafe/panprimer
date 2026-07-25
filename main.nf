#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Pangenome PCR primer design — two-stage pipeline with cache-aware projection.
//   ANCHOR -> PROJECT (per hap, cache-aware) -> MASK_DESIGN
//          -> STAGE_A (coverage shortlist) -> EVALUATE (per hap, genome-wide on shortlist)
//          -> AGGREGATE
// Each process shells out to a `pangenome-primer` subcommand. Prebuilt per-haplotype caches
// (.fai/.gzi for BGZF random access, .anchors.tsv.gz for projection) are staged alongside
// the FASTA so nothing is rebuilt inside a work dir. Whatever is absent is simply not
// staged, so a lean checkout stages ~0.91 GB per haplotype instead of ~15 GB.

process ANCHOR {
    input:
      val target
      path chm13
      path cfg
      val assembly
      path grch38          // assets/NO_FILE placeholder when not using GRCh38 coords
    output:
      path 'anchor.json'
    script:
    def g = grch38.name != 'NO_FILE' ? "--grch38 ${grch38}" : ''
    """
    pangenome-primer anchor --target '${target}' --chm13 ${chm13} \
        --target-assembly ${assembly} ${g} --config ${cfg} --out anchor.json
    """
}

process PROJECT {
    tag "${hid}"
    input:
      path anchor
      tuple val(hid), path(fa), path(idx)
    output:
      tuple val(hid), path("proj.${hid.replace('#','_')}.json")
    script:
    """
    pangenome-primer project --anchor ${anchor} --hap-fasta ${fa} \
        --hap-id '${hid}' --out proj.${hid.replace('#','_')}.json
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

process STAGE_A {
    input:
      path candidates
      path projections
      val mode
      val top_k
    output:
      path 'shortlist.json'
    script:
    def proj_args = projections.collect { "--projection ${it}" }.join(' ')
    def mode_arg = mode ? "--mode ${mode}" : ''
    """
    pangenome-primer stage-a --candidates ${candidates} ${proj_args} \
        ${mode_arg} --top-k ${top_k} --out shortlist.json
    """
}

process EVALUATE {
    tag "${hid}"
    input:
      path shortlist
      tuple val(hid), path(proj), path(fa), path(idx)
      val mode
    output:
      path "result.${hid.replace('#','_')}.json"
    script:
    def mode_arg = mode ? "--mode ${mode}" : ''
    // This is the fan-out: one task per haplotype, many at once. Left alone, the compiled
    // scanner sizes its pool from the HOST's core count, so on a 256-core box each of 64
    // concurrent tasks spawns 256 workers -- 16,384 threads over 256 cores. Handing it
    // `task.cpus` makes threads x concurrency <= the allocation by construction rather than
    // by convention. Tune the number via `withName: EVALUATE { cpus = N }` in
    // nextflow.config; see ISSUE-002 in docs/issues.md for the measured curve.
    """
    export PGP_SCAN_THREADS=${task.cpus}
    pangenome-primer evaluate --candidates ${shortlist} --projection ${proj} \
        --hap-fasta ${fa} ${mode_arg} --out result.${hid.replace('#','_')}.json
    """
}

process AGGREGATE {
    publishDir "${params.outdir}", mode: 'copy'
    input:
      path anchor
      path shortlist
      path results
    output:
      path 'report/*'
    script:
    def res_args = results.collect { "--result ${it}" }.join(' ')
    """
    pangenome-primer aggregate --anchor ${anchor} --candidates ${shortlist} \
        ${res_args} --outdir report
    """
}

workflow {
    if( !params.target || !params.chm13 )
        error "Required: --target 'chr:start-end'|FASTA and --chm13 <CHM13 FASTA>"

    cfg    = file(params.config)
    chm13  = file(params.chm13)
    mode   = params.mode ?: ''
    top_k  = params.top_k ?: 5
    grch38 = file(params.grch38 ?: "${projectDir}/assets/NO_FILE")

    // haplotype channel: (hap_id, fasta, [prebuilt caches]) from samples.tsv
    haplos = Channel.fromPath(params.samples)
        .splitCsv(sep: '\t', header: false)
        .filter { it && !it[0].startsWith('#') && it[0] != 'sample' }
        .map { row ->
            def fa = file(row[3])
            // Sidecars may be named for the .fa.gz itself or for the .gz-stripped stem
            // -- collect both and dedupe, mirroring `samples.sidecar_path`.
            // `exists()` is load-bearing: for a pattern with no glob metacharacter `files()`
            // hands back the literal path whether or not anything is there, so without this
            // filter every optional sidecar (e.g. a PAF that was never built) would be
            // staged as a broken symlink and fail the task.
            def stem = fa.toString().replaceFirst(/\.gz$/, '')
            def idx = ['fai','gzi','anchors.tsv.gz','mmi','chm13.paf']
                        .collectMany { files("${fa}.${it}") + files("${stem}.${it}") }
                        .findAll { it.exists() }
                        .unique()
            tuple("${row[0]}#hap${row[1]}", fa, idx)
        }

    anchor    = ANCHOR(params.target, chm13, cfg, params.target_assembly ?: 'chm13', grch38)
    projected = PROJECT(anchor, haplos)                     // (hid, proj.json)
    projjson  = projected.map { it[1] }.collect()
    candidates = MASK_DESIGN(anchor, projjson)
    shortlist = STAGE_A(candidates, projjson, mode, top_k)  // top-K by coverage
    evalIn    = projected.join(haplos)                      // (hid, proj.json, fa, idx)
    results   = EVALUATE(shortlist, evalIn, mode)           // genome-wide on shortlist only
    AGGREGATE(anchor, shortlist, results.collect())
}
