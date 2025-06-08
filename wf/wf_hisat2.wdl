version 1.0

workflow wf_hisat2 {
  input {
    File read1
    File? read2
    String samplename
    File? index
    File? genome_fasta
  }

  if (!defined(index)) {
    call hisat2_index {
      input:
        genome_fasta = genome_fasta
    }
  }

  if (defined(index)) {
    call extract_index {
      input:
        salmon_index_tar = index
    }
  }

  call hisat2_align {
    input:
      hisat2_index_files = select_first([extract_index.index_files, hisat2_index.index_files]),
      read1 = read1,
      read2 = read2,
      samplename = samplename
  }

  output {
    File aligned_sam = hisat2_align.sam
  }
}


task hisat2_index {
  input {
    File genome_fasta
  }

  command <<< 
    mkdir -p hisat2_index
    hisat2-build ~{genome_fasta} hisat2_index/genome
  >>>

  output {
    Array[File] index_files = glob("hisat2_index/*")
  }

  runtime {
    docker: "quay.io/biocontainers/hisat2:2.2.1--py38h7f98852_0"
    cpu: 4
    memory: "8G"
  }
}

task hisat2_align {
  input {
    Array[File] hisat2_index_files
    File read1
    File? read2
    String samplename
  }

  command <<< 
    mkdir -p hisat2_index
    for f in ~{sep=' ' hisat2_index_files}; do
      cp "$f" hisat2_index/
    done

    if [ -n "~{read2}" ]; then
      hisat2 -x hisat2_index/genome -1 ~{read1} -2 ~{read2} -S ~{samplename}.sam
    else
      hisat2 -x hisat2_index/genome -U ~{read1} -S ~{samplename}.sam
    fi
  >>>

  output {
    File sam = "~{samplename}.sam"
  }

  runtime {
    docker: "quay.io/biocontainers/hisat2:2.2.1--py38h7f98852_0"
    cpu: 4
    memory: "8G"
  }
}
