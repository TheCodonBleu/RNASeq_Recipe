version 1.0

task feature_counts {
  input {
    File? bam
    File? gtf
  }

  command <<<
    featureCounts -a ~{gtf} -o counts.txt ~{bam}
  >>>

  output {
    File? counts = "counts.txt"
  }

  runtime {
    docker: "biocontainers/subread:2.0.1--h7132678_0"
  }
}
