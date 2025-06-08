version 1.0

task decompress_file {
  input {
    File? file
  }

  command <<<
    set -e

    filename=$(basename "~{file}")

    # Extract file extension (handling multiple extensions like .gtf.gz)
    if [[ "$filename" == *.gz ]]; then
      echo "Detected .gz compressed file"
      orig_name=$(echo "$filename" | sed 's/\.gz$//')
      gunzip -c "~{file}" > "$orig_name"

    elif [[ "$filename" == *.zip ]]; then
      echo "Detected .zip compressed file"
      # Extract the first file (assumes zip has only one)
      unzip -p "~{file}" > temp_out
      orig_name=$(unzip -l "~{file}" | awk 'NR==4{print $4}')
      mv temp_out "$orig_name"

    else
      echo "Detected uncompressed file"
      cp "~{file}" "$filename"
      orig_name="$filename"
    fi

    echo "Decompressed to: $orig_name"
  >>>

  output {
    File file_decompressed = glob("*")[0]
  }

  runtime {
    docker: "ubuntu:20.04"
  }
}
