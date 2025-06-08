version 1.0

import '../tasks/utils/decompress.wdl' as Decompress

workflow wf_salmon {
  input {
    File read1
    File read2  
    String samplename
    File? index
    File? transcripts
  }

  # If no index provided, create the index using the transcript fasta file (compressed)
  if (!defined(index)) {

    call salmon_index {
      input:
        transcripts = transcripts
    } 
  }

  if (defined(index)) {
    call extract_index {
      input:
        salmon_index_tar = index
    }
  }

  call salmon_quant {
    input:
      salmon_index_files = select_first([extract_index.index_files, salmon_index.index_files]),
      read1 = read1,
      read2 = read2,
      samplename = samplename
  }

  output {
    File salmon_quant_sf = salmon_quant.quant_sf
    File salmon_lib_format_counts = salmon_quant.lib_format_counts
    File salmon_cmd_info = salmon_quant.cmd_info
    File salmon_quant_log = salmon_quant.log

  }
}

task extract_index {
  input {
    File? salmon_index_tar
  }

  command <<<
    mkdir -p salmon_index
    echo "Index found, dezip :"
    echo "execute : tar -xvf ~{salmon_index_tar} -C ."
    tar -xvf ~{salmon_index_tar} -C salmon_index
  >>>

  output {
    Array[File] index_files = glob("salmon_index/*/*")
  }

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "2G"
  }
}

task salmon_quant {
  input {
    Array[File] salmon_index_files
    File read1
    File read2  
    String samplename
  }

  command <<<
    mkdir index
    for f in ~{sep=' ' salmon_index_files}; do
      cp "$f" index/
    done
    
    if [[ -n "~{read2}" ]]; then
      echo "Paired-end mode"
      salmon quant -i index -l A \
        -1 ~{read1} \
        -2 ~{read2} \
        -p 8 --validateMappings -o salmon_out
    else
      echo "Single-end mode"
      salmon quant -i index -l A \
        -1 ~{read1} \
        -p 8 --validateMappings -o salmon_out

    fi

    cp salmon_out/quant.sf ~{samplename}.quant.sf
    cp salmon_out/lib_format_counts.json ~{samplename}.lib_format_counts.json
    cp salmon_out/cmd_info.json ~{samplename}.cmd_info.json
    cp salmon_out/logs/salmon_quant.log ~{samplename}.salmon_quant.log

  >>>

  output {
    File quant_sf = "~{samplename}.quant.sf"
    File lib_format_counts = "~{samplename}.lib_format_counts.json"
    File cmd_info = "~{samplename}.cmd_info.json"
    File log = "~{samplename}.salmon_quant.log"
  }

  runtime {
    docker: "quay.io/biocontainers/salmon:1.9.0--h7e5ed60_0"
    cpu: 8
    memory: "8G"
  }
}

task salmon_index {
  input {
    File? transcripts
  }

  command <<<
    mkdir -p salmon_index

    salmon index \
      -t ~{transcripts} \
      -i salmon_index \
  >>>

  output {
    Array[File] index_files = glob("salmon_index/*")
  }

  runtime {
    docker: "combinelab/salmon:1.10.2"
    cpu: 25
    memory: "4G"
  }
}
