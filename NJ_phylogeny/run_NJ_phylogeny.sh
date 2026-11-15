#!/bin/bash
#SBATCH --job-name=timopheevii_NJ
#SBATCH --output=NJ_phylogeny_%j.out
#SBATCH --error=NJ_phylogeny_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=800G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

# Rooted bootstrapped NJ tree of T. timopheevii (KAUST IBEX).
# Step 1 (VCF -> GDS) runs only if the GDS file does not exist yet.
# Usage: sbatch run_NJ_phylogeny.sh

module purge
module load R

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

if [[ ! -f timopheevii_WGS_SNPs.gds ]]; then
    Rscript 01_vcf_to_gds.R
fi

Rscript 02_bootstrap_NJ_tree.R
