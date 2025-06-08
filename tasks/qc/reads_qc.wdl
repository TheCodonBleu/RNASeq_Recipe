version 1.0

task fastqc {
  input {
    File read1
    File read2

    Int memory = 4
    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/fastqc:0.12.1"
  }
  String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
  
  command <<<
      # get fastqc version
    fastqc --version | tee VERSION
    
  
    fastqc \
      --outdir . \
      --threads ~{cpu} \
      --extract \
      ~{read1} \
      ~{read2}

    grep "Total Sequences" ~{read1_name}_fastqc/fastqc_data.txt | cut -f 2 | tee READ1_SEQS
    read1_seqs=$(cat READ1_SEQS)
    grep "Total Sequences" ~{read2_name}_fastqc/fastqc_data.txt | cut -f 2 | tee READ2_SEQS
    read2_seqs=$(cat READ2_SEQS)

    # capture number of read pairs
    if [ "${read1_seqs}" == "${read2_seqs}" ]; then
      read_pairs=${read1_seqs}
    else
      read_pairs="Uneven pairs: R1=${read1_seqs}, R2=${read2_seqs}"
    fi
    
    echo "$read_pairs" | tee READ_PAIRS
  >>>

  output {
    File read1_fastqc_html = "~{read1_name}_fastqc.html"
    File read1_fastqc_zip = "~{read1_name}_fastqc.zip"
    File read2_fastqc_html = "~{read2_name}_fastqc.html"
    File read2_fastqc_zip = "~{read2_name}_fastqc.zip"
    
    String fastqc_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}

task multiqc {
  input {
    Array[File] input_dirs
  }

  command {
    multiqc ~{sep=' ' input_dirs}
  }

  output {
    File report = "multiqc_report.html"
  }

  runtime {
    docker: "ewels/multiqc:1.11"
  }
}


