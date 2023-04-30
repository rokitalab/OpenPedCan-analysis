cwlVersion: v1.2
class: CommandLineTool
id: liftover-collapse-rnaseq
doc: "Tool to use a GTF and expression matrix with ENSG to liftover gene symbols and collapse repeat symbols if needed"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'python:3.9'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: liftover_collapse_rnaseq.py
        entry:
          $include: ../scripts/liftover_collapse_rnaseq.py


baseCommand: [python3 liftover_collapse_rnaseq.py]

inputs:
  table: {type: File, doc: "Input expression matrix. Can be gzipped or plain tsv",
    inputBinding: { prefix: "--table", position: 1 } }
  output_basename: {type: string, doc: "output basename of liftover and collapse if flag given",
    inputBinding: { prefix: "--output-basename", position: 1 } }
  input_gtf: {type: File, doc: "GTF file to use as liftover reference. Can be gzipped or plain tsv",
    inputBinding: { prefix: "--input-gtf", position: 1 } }
  gene_id: {type: string, doc: "NAME of gene ID (ENSG values) column",
    inputBinding: { prefix: "--gene-id", position: 1 } }
  gene_name: {type: string, doc: "NAME of gene name column, if present to keep if not found",
    inputBinding: { prefix: "--gene-name", position: 1 } }
  skip: {type: 'int?', doc: "Number of lines to skip if needed, i.e. file has extra unneeded header lines",
    inputBinding: { prefix: "--skip", position: 1 } }
  collapse: {type: 'boolean?', doc: "If set, will collapse on repeat gene symbols and choose highest mean expression", default: true,
    inputBinding: { prefix: "--collapse", position: 1 } }




outputs:
  liftover:
    type: File
    outputBinding:
      glob: '*.liftover.tsv'
    doc: "Un-collapsed expression matrix with lifted-over gene symbols"
  collapsed:
    type: File
    outputBinding:
      glob: '*.collapsed.tsv'
    doc: "Collapsed expression matrix with lifted-over gene symbols"
 