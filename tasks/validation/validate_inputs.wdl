version 1.0

task validate_aligner {
  
  input {
    String aligner
  }

  command <<<
    if [ "~{aligner}" != "salmon" ] && [ "~{aligner}" != "star" ]; then
      echo "ERROR: aligner must be 'salmon' or 'star'" >&2
      exit 1
    fi
  >>>


  runtime {
    docker: "ubuntu:20.04"
  }
}