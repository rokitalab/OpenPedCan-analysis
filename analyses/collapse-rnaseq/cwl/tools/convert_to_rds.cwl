cwlVersion: v1.2
class: CommandLineTool
id: convert-to-rds
doc: "Converts tsv table to rds"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: convert_to_rds.R
        entry:
          $include: ../scripts/convert_to_rds.R
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)

baseCommand: [Rscript convert_to_rds.R]

inputs:
  input_tsv: {type: File, doc: "Input expression matrix. Can be gzipped or plain tsv",
    inputBinding: { prefix: "--input-tsv", position: 1} }
  gene_col: {type: 'string?', doc: "Name of column with gene symbol", default: 'gene_name',
    inputBinding: { prefix: "--gene-col", position: 1} }
  output_filename: {type: string, doc: "Output file name, should fit pattern gene-[expression/counts]-rsem-[fpkm/tpm/expected_count]-collapsed.rds",
    inputBinding: { prefix: "--output-filename", position: 1} }
  ram: { type: 'int?', doc: "Set ram requirement. Unfortunately reading and writing rds files is beefy!", default: 32 }

outputs:
  pirate_output:
    type: File
    outputBinding:
      glob: '*.rds'
    doc: "rds file generated from input tsv"
