cwlVersion: v1.2
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
    - $(inputs.input_files)

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      echo '${return inputs.input_files.map(function(e){return e.path}).join('\n')}' > files.txt
      && cat files.txt | xargs -IFN -P $(inputs.cores) gzip FN

inputs:
  input_files: { type: 'File[]', doc: "Files to gzip" }
  cores: { type: 'int?', doc: "Num threads to use. Uses xargs to gzip files in parallel", default: 2 }

outputs:
  gzipped_files:
    type: 'File[]'
    outputBinding:
      glob: '*.gz'
    doc: "Gzipped inputs"
  file_list:
    type: File
    outputBinding:
      glob: 'files.txt'
    doc: "More of a debug file"
