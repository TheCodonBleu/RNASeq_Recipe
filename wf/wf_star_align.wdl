version 1.0

workflow star_align_wf {
  input {
    File read1
    File read2
    File? index_compressed
    File? ref_genome_fasta
    File? gtf_file
  }

  if (defined(index_compressed)) {
    call star_extract_index {
      input:
        index_compressed = index_compressed
    }

    call star_align as star_align_index_compressed {
      input:
        read1 = read1,
        read2 = read2,
        index_files = star_extract_index.index_files
    }
  }

  if (!defined(index_compressed)) {
    call star_index {
      input:
        ref_genome_fasta = ref_genome_fasta,
        gtf_file = gtf_file
    }

    call star_align as star_align_index_generated {
      input:
        read1 = read1,
        read2 = read2,
        index_files = star_index.index_files
    }
  }

  output {
    File bam = select_first([star_align_index_compressed.bam, star_align_index_generated.bam])
  }
}

task star_align {
  input {
    File read1
    File read2
    Array[File] index_files
  }

    command <<<
    mkdir star_index
    for f in ~{sep=' ' index_files}; do
        cp "$f" star_index/
    done

    STAR \
        --genomeDir star_index \
        --readFilesIn ~{read1} ~{read2} \
        --runThreadN 4 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix star_ 2> star_error.log
    >>>

  output {
    File? bam = "star_Aligned.sortedByCoord.out.bam"
  }

  runtime {
    docker: "alexdobin/star:2.7.10a_alpha_220506"
    cpu: 4
    memory: "8G"
  }
}

task star_extract_index {
  input {
    File? index_compressed
  }

  command {
    mkdir star_index
    if [[ "~{index_compressed}" == *.zip ]]; then
      unzip ~{index_compressed} -d star_index
    else
      tar -xvf ~{index_compressed} -C star_index
    fi
  }

  output {
    Array[File] index_files = glob("star_index/*")
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "4G"
    cpu: 1
  }
}

task star_index {
  input {
    File? ref_genome_fasta
    File? gtf_file

    Int cpu = 25
  }

  command <<<
    mkdir star_index

    if [[ "~{ref_genome_fasta}" == *.gz ]]; then
      echo "Decompressing FASTA file..."
      gunzip -c ~{ref_genome_fasta} > genome.fa
      fasta_file="genome.fa"
    else
      fasta_file="~{ref_genome_fasta}"
    fi
    
    STAR \
      --runMode genomeGenerate \
      --runThreadN ~{cpu} \
      --genomeDir star_index \
      --genomeFastaFiles ${fasta_file} \
      --sjdbGTFfile ~{gtf_file} \
      --sjdbOverhang 100
  >>>

  output {
    Array[File] index_files = glob("star_index/*")
  }

  runtime {
    docker: "alexdobin/star:2.7.10a_alpha_220506"
    memory: "8G"
    cpu: cpu
  }
}
