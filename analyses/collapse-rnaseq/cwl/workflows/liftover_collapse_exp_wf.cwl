cwlVersion: v1.2
class: Workflow
id: liftover-collapse-rnaseq-wf
label: Liftover and Collapse RNAseq WF
doc: "Lifts over gene symbols from given gtf and output lifted and collapsed by gene symbol tsv, and rds of collapsed file"


inputs:
  # lift and collapse
  table: { type: File, doc: "Input expression matrix. Can be gzipped or plain tsv" }
  output_basename:  {type: string, doc: "output basename of liftover and collapse if flag given" }
  input_gtf: { type: File, doc: "GTF file to use as liftover reference. Can be gzipped or plain tsv" }
  gene_id: { type: string, doc: "NAME of gene ID (ENSG values) column"  }
  gene_name: { type: string, doc: "NAME of gene name column, if present to keep if not found" }
  skip: { type: 'int?', doc: "Number of lines to skip if needed, i.e. file has extra unneeded header lines" }
  collapse: { type: 'boolean?', doc: "If set, will collapse on repeat gene symbols and choose highest mean expression", default: true }
  # convert to rds
  output_filename: { type: string, doc: "Output file name, should fit pattern gene-[expression/counts]-rsem-[fpkm/tpm/expected_count]-collapsed.rds" }
  rds_ram: { type: 'int?', doc: "Set ram requirement. Unfortunately reading and writing rds files is beefy!", default: 32}

outputs:
  gzipped_results: { type: 'File[]', outputSource: gzip_tsv/gzipped_files }
  expression_rds: { type: 'File?', outputSource: convert_to_rds/pirate_output }

steps:
  liftover_collapse_rnaseq:
    run: ../tools/liftover_collapse_rnaseq.cwl
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.2xlarge;ebs-gp2;400
      doc: "Chosen for speed and lower cost"
    in:
      table: table
      output_basename: output_basename
      input_gtf: input_gtf
      gene_id: gene_id
      gene_name: gene_name
      skip: skip
      collapse: collapse
    out: [liftover, collapsed]

  gzip_tsv:
    run: ../tools/ubuntu_gzip.cwl
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.2xlarge;ebs-gp2;400
      doc: "Chosen for speed and lower cost"
    in:
      input_files:
        source: [liftover_collapse_rnaseq/liftover, liftover_collapse_rnaseq/collapsed]
        valueFrom: |
          $([self])
    out: [gzipped_files]

  convert_to_rds:
    run: ../tools/convert_to_rds.cwl
    when: $(inputs.input_tsv != null)
    in:
      input_tsv: liftover_collapse_rnaseq/collapsed
      output_filename: output_filename
      ram: rds_ram
    out: [pirate_output]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
