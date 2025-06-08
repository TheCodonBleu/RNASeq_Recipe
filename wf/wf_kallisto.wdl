version 1.0

workflow wf_kallisto {
  input {
    File read1
    File read2
    String samplename
    File? index
    File? transcripts 
  }


  if (!defined(index)) {
    call kallisto_index {
      input:
        transcripts = transcripts
    }
  }


  call kallisto_quant {
    input:
      kallisto_index_files = select_first([kallisto_index.index_files, []]),
      read1 = read1,
      read2 = read2,
      samplename = samplename
  }

  output {
    File kallisto_quant_abundance = kallisto_quant.abundance
    File kallisto_run_info = kallisto_quant.run_info
  }
}

task kallisto_index {
  input {
    File? transcripts
  }

  command <<< 
    mkdir -p kallisto_index
    kallisto index -i kallisto_index/transcripts.idx ~{transcripts}
  >>>

  output {
    Array[File] index_files = glob("kallisto_index/*")
  }

  runtime {
    docker: "danhumassmed/salmon-kallisto:1.0.1"
    cpu: 4
    memory: "4G"
  }
}

task kallisto_quant {
  input {
    Array[File] kallisto_index_files
    File read1
    File read2
    String samplename
  }

  command <<< 
    mkdir index
    for f in ~{sep=' ' kallisto_index_files}; do
      cp "$f" index/
    done

    echo "Kallisto Quantification"

    kallisto quant -i index/transcripts.idx -o kallisto_out -t 8 -b 100 ~{read1} ~{read2}
    
    cp kallisto_out/abundance.tsv ~{samplename}.abundance.tsv
    cp kallisto_out/run_info.json ~{samplename}.run_info.json
  >>>

  output {
    File abundance = "~{samplename}.abundance.tsv"
    File run_info = "~{samplename}.run_info.json"
  }

  runtime {
    docker: "danhumassmed/salmon-kallisto:1.0.1"
    cpu: 8
    memory: "8G"
  }
}
