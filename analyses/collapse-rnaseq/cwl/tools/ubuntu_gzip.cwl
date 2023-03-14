cwlVersion: v1.0
class: CommandLineTool
id: ubuntu-gzip
doc: "Simple tool to gzip files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:22.04'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entryname: files.txt
      entry: |
        $(inputs.input_files.map(function(e){return e.path}).join('\n'))

baseCommand: [cat]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      files.txt | xargs -IFN -P $(inputs.cores) gzip FN

inputs:
  input_files: { type: 'File[]', doc: "Files to gzip" }
  cores: { type: 'int?', doc: "Num threads to use. Uses xargs to gzip files in parallel", default: 2 }

outputs:
  gzipped_files:
    type: 'File[]'
    outputBinding:
      glob: '*.gz'
    doc: "Gzipped inputs"
