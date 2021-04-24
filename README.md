# Data preprocessing for variant discovery

Data preprocessing workflow for variant discovery on UT HPC. Based on GATK Best Practice: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery.

## Dependencies

* Cromwell
* Conda:
    * BWA
    * Samtools
    * GATK4

Environment: [preprocessing.yaml](preprocessing.yaml)

## Usage

1. Create `inputs.json` file for preprocessing using:
```bash
./input.py [-h] path
```
2. Run preprocessing in SLURM using:
```bash
sbatch cromwell.sh
```
