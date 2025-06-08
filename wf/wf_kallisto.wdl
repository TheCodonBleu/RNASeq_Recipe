version 1.0

workflow wf_kallisto {
  input {
    File read1
    File read2
    String samplename
    File? kallisto_index
    File? transcripts 
  }


  if (!defined(kallisto_index)) {
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
    File abundance_tsv = kallisto_quant.abundance
  }
}

task kallisto_index {
  input {
    File transcripts
  }

  command <<< 
    mkdir -p kallisto_index
    kallisto index -i kallisto_index/transcripts.idx ~{transcripts}
  >>>

  output {
    Array[File] index_files = glob("kallisto_index/*")
  }

  runtime {
    docker: "biocontainers/kallisto:v0.46.2_cv1"
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

    kallisto quant -i index/transcripts.idx -o kallisto_out -t 8 -b 100 -1 ~{read1} -2 ~{read2}
    
    cp kallisto_out/abundance.tsv ~{samplename}.kallisto.abundance.tsv
  >>>

  output {
    File abundance = "~{samplename}.kallisto.abundance.tsv"
  }

  runtime {
    docker: "biocontainers/kallisto:v0.46.2_cv1"
    cpu: 8
    memory: "8G"
  }
}
