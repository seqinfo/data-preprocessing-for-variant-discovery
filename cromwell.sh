#!/bin/bash

#SBATCH --job-name=preprocessing
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/cromwell/%j.log

module load java-1.8.0_40


java -Dconfig.file=slurm.conf -Xmx5g -jar cromwell-59.jar run -i inputs.json \
  -o options.json workflows/main.wdl
