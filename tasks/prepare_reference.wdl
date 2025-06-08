version 1.0

task PrepareReference {
  input {
    File fasta
    File gtf
    String species
  }

  command <<<
    echo "Preparing reference for ~{species}..."
    mkdir ref
    cp ~{fasta} ref/genome.fa
    cp ~{gtf} ref/genes.gtf
    salmon index -t ref/genome.fa -g ref/genes.gtf -i ref/salmon_index
  >>>

  output {
    File ref_dir = "ref"
  }

  runtime {
    docker: "compsci/salmon:latest"
    memory: "8G"
    cpu: 2
  }
}
