version 1.0

import "preprocessing.wdl" as Preprocessing

workflow Main {
  input {
    String library_id
    Array[String] sample_ids
    Map[String, Array[File]] samples
  }

  scatter (sample_id in sample_ids) {
    Array[File] sample_lanes = samples[sample_id]

    call Preprocessing.Preprocessing as Preprocessing {
      input:
        library_id = library_id,
        sample_id = sample_id,
        sample_lanes = sample_lanes
    }
  }

  output {
    Array[File] bams = Preprocessing.bqsr_bam
    Array[File] bams_idx = Preprocessing.bqsr_bam_idx
  }
}
