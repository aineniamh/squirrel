process squirrel {

  container "${params.wf.container}:${workflow.manifest.version}"
  
  publishDir "${params.out_dir}", mode: 'copy', saveAs: {  fn -> fn.replace("squirrel_output/", "")}

  input:
    path fasta
    path refs
    path additional_mask
    
  output:
    path "squirrel_output/${fasta.baseName}.aln.fasta"
    path "squirrel_output/**"

  script:
    extra = ""
    if ( params.output_intermediates )
        extra += " --no-temp"
    if ( params.seq_qc )
        extra += " --seq-qc"
    if ( params.assembly_refs )
        extra += " --assembly-refs ${refs}"
    if ( params.no_mask )
        extra += " --no-mask"
    if ( params.no_itr_mask )
        extra += " --no-itr-mask"
    if ( params.additional_mask )
        extra += " --additional-mask ${additional_mask}"
    if ( params.extract_cds )
        extra += " --extract-cds"
    if ( params.concatenate )
        extra += " --concatenate"
    if ( params.clade )
        extra += " --clade ${params.clade}"
    if ( params.run_phylo )
        extra += " --run-phylo"
    if ( params.outgroups )
        extra += " --outgroups ${params.outgroups}"
    
    """
    squirrel ${fasta} -o squirrel_output --outfile ${fasta.baseName}.aln.fasta --tempdir squirrel_tmp -t ${task.cpus} ${extra}
    """
}

workflow {
  if ( params.assembly_refs ) {
      refs_ch = channel.fromPath("${params.assembly_refs}", checkIfExists:true)
  } else {
      refs_ch = channel.fromPath("${projectDir}/${params.default_assembly_refs}", checkIfExists:true)
  }

if ( params.additional_mask ) {
      mask_ch = channel.fromPath("${params.additional_mask}", checkIfExists:true)
  } else {
      mask_ch = channel.fromPath("${projectDir}/${params.default_additional_mask}", checkIfExists:true)
  }
  
  fasta_ch = Channel.of(file("${params.fasta}", type: "file", checkIfExists:true))
  
  squirrel(fasta_ch, refs_ch, mask_ch)
}
