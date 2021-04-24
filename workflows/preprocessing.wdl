version 1.0

workflow Preprocessing {
  input {
    String library_id
    String sample_id
    Array[File] sample_lanes

    File ref_fasta = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta"
    File ref_fasta_idx = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.fai"
    File ref_dict = "gatk_resource_bundle/Homo_sapiens_assembly38.dict"
    File ref_alt = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.alt"
    File ref_amb = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.amb"
    File ref_ann = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.ann"
    File ref_bwt = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.bwt"
    File ref_pac = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.pac"
    File ref_sa = "gatk_resource_bundle/Homo_sapiens_assembly38.fasta.64.sa"

    Array[File] known_sites_vcfs = [
      "dbsnp154/dbsnp_154.hg38.vcf.gz"
    ]
    Array[File] known_sites_vcf_idxs = [
      "dbsnp154/dbsnp_154.hg38.vcf.gz.tbi"
    ]

    String conda_env
  }

  scatter (sample_lane in sample_lanes) {
    String lane_id = basename(sample_lane, ".fastq.gz")

    call AlignSort {
      input:
        library_id = library_id,
        sample_id = sample_id,
        lane_id = lane_id,
        fastq = sample_lane,
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        conda_env = conda_env
    }
  }

  call MarkDuplicatesSort {
    input:
      sample_id = sample_id,
      bams = AlignSort.bam,
      conda_env = conda_env
  }

  call FixTags {
    input:
      sample_id = sample_id,
      bam = MarkDuplicatesSort.markdups_bam,
      bam_idx = MarkDuplicatesSort.markdups_bam_idx,
      ref_fasta = ref_fasta,
      ref_fasta_idx = ref_fasta_idx,
      ref_dict = ref_dict,
      conda_env = conda_env
  }

  call BaseRecalibratorSpark {
    input:
      sample_id = sample_id,
      bam = FixTags.fixtags_bam,
      bam_idx = FixTags.fixtags_bam_idx,
      ref_fasta = ref_fasta,
      ref_fasta_idx = ref_fasta_idx,
      ref_dict = ref_dict,
      known_sites_vcfs = known_sites_vcfs,
      known_sites_vcf_idxs = known_sites_vcf_idxs,
      conda_env = conda_env
  }

  call ApplyBQSR {
    input:
      sample_id = sample_id,
      bam = FixTags.fixtags_bam,
      bam_idx = FixTags.fixtags_bam_idx,
      recal_table = BaseRecalibratorSpark.recal_table,
      ref_fasta = ref_fasta,
      ref_fasta_idx = ref_fasta_idx,
      ref_dict = ref_dict,
      conda_env = conda_env
  }

  output {
    File bqsr_bam = ApplyBQSR.bqsr_bam
    File bqsr_bam_idx = ApplyBQSR.bqsr_bam_idx
  }
}

task AlignSort {
  input {
    String library_id
    String sample_id
    String lane_id
    File fastq
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String conda_env
  }

  Int cpus = 10

  command <<<
    source activate ~{conda_env}
    set -euxo pipefail

    # make read groups
    header="$(zgrep -m 1 "^@" ~{fastq})"
    flowcell="$(echo "${header}" | cut -d ":" -f 3)"
    lane="$(echo "${header}" | cut -d ":" -f 4)"
    index="$(echo "${header}" | cut -d ":" -f 10)"
    id="${flowcell}"."${lane}"
    pu="${flowcell}"."${lane}"."${index}"
    sm=~{sample_id}
    pl=ILLUMINA
    lb=~{library_id}

    # align and sort
    bwa mem -t ~{cpus} -K 100000000 \
      -R "@RG\tID:"${id}"\tPU:"${pu}"\tSM:"${sm}"\tPL:"${pl}"\tLB:"${lb}"" \
      ~{ref_fasta} ~{fastq} \
      | samtools sort -n -o ~{lane_id}_qsort.bam -T /tmp -@ ~{cpus} -
  >>>

  output {
    File bam = "~{lane_id}_qsort.bam"
  }

  runtime {
    time: 30
    cpus: cpus
    memory: 12
  }
}

task MarkDuplicatesSort {
  input {
    String sample_id
    Array[File] bams
    String conda_env
  }

  command {
    source activate ~{conda_env}

    gatk --java-options "-Xmx36g" MarkDuplicatesSpark -I ~{sep=" -I " bams} \
      -O ~{sample_id}_markdups.bam --tmp-dir /tmp
  }

  output {
    File markdups_bam = "~{sample_id}_markdups.bam"
    File markdups_bam_idx = "~{sample_id}_markdups.bam.bai"
  }

  runtime {
    time: 30
    cpus: 10
    memory: 36
  }
}

task FixTags {
  input {
    String sample_id
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    String conda_env
  }

  command {
    source activate ~{conda_env}

    gatk --java-options "-Xmx11g" SetNmMdAndUqTags -I ~{bam} \
      -O ~{sample_id}_fixtags.bam -R ~{ref_fasta} --CREATE_INDEX --TMP_DIR /tmp
  }

  runtime {
    time: 30
    cpus: 1
    memory: 12
  }

  output {
    File fixtags_bam = "~{sample_id}_fixtags.bam"
    File fixtags_bam_idx = "~{sample_id}_fixtags.bai"
  }
}

task BaseRecalibratorSpark {
  input {
    String sample_id
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    Array[File] known_sites_vcfs
    Array[File] known_sites_vcf_idxs
    String conda_env
  }

  command {
    source activate ~{conda_env}

    gatk --java-options "-Xmx41g" BaseRecalibratorSpark -I ~{bam} \
      --known-sites ~{sep=" --known-sites " known_sites_vcfs} \
      -O ~{sample_id}_recal.table -R ~{ref_fasta} --tmp-dir /tmp
  }

  runtime {
    time: 60
    cpus: 10
    memory: 42
  }

  output {
    File recal_table = "~{sample_id}_recal.table"
  }
}

task ApplyBQSR {
  input {
    String sample_id
    File bam
    File bam_idx
    File recal_table
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    String conda_env
  }


  command {
    source activate ~{conda_env}

    gatk --java-options "-Xmx5g" ApplyBQSR -bqsr ~{recal_table} -I ~{bam} \
      -O ~{sample_id}_bqsr.bam -R ~{ref_fasta} --tmp-dir /tmp
  }

  runtime {
    time: 30
    cpus: 1
    memory: 6
  }

  output {
    File bqsr_bam = "~{sample_id}_bqsr.bam"
    File bqsr_bam_idx = "~{sample_id}_bqsr.bai"
  }
}
