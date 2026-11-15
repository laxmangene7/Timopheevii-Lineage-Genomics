# Rooted bootstrapped NJ phylogeny of *T. timopheevii*

Builds a rooted neighbour-joining (NJ) tree with 1,000 bootstrap replicates from the *T. timopheevii* WGS SNP set. The tree is rooted on wild emmer (*T. dicoccoides*).

## Scripts

| Script | What it does |
| --- | --- |
| `run_NJ_phylogeny.sh` | SLURM job for KAUST IBEX (10 CPUs, 800 GB RAM, 24 h, `batch` partition). It runs step 1 only if the GDS file doesn't exist yet, then runs step 2. |
| `01_vcf_to_gds.R` | Converts `timopheevii_WGS_SNPs.vcf.gz` to GDS format (biallelic SNPs). Prints sample and SNP counts and checks that the outgroup accessions are present. |
| `02_bootstrap_NJ_tree.R` | Computes pairwise distances, the reference NJ tree, 1,000 bootstrap trees and node support, then roots the tree on the outgroup. |

## Requirements

R with `ape` (CRAN) and `SNPRelate` and `gdsfmt` (Bioconductor).

## Input

`timopheevii_WGS_SNPs.vcf.gz` (with its `.csi` index). The VCF is too large for GitHub and is available on Dryad: [doi:10.5061/dryad.ksn02v7g7](https://doi.org/10.5061/dryad.ksn02v7g7). Download it into this folder before running. This is the SNP set from reads mapped to the *T. timopheevii* TA2804 v1.3 reference. Steps that produced it, as recorded in the VCF header:

1. `bcftools mpileup` (SNPs only) and `bcftools call -m`
2. `bcftools filter`: QUAL ≥ 30, INFO/DP ≥ 20, ≤ 20% missing, 0.01 < AF < 0.99
3. Random 10% of SNPs per chromosome (GATK `SelectVariants --select-random-fraction 0.1`)
4. Chromosomes concatenated with `bcftools concat`

## Method

- **Samples:** IG46434, IG116164 and IG116166 were excluded because we identified them as non-*timopheevii* accessions. SNPs that are monomorphic among the remaining samples are dropped.
- **Distance:** For each sample pair, allele sharing is `1 - |g_i - g_j| / 2`, averaged over SNPs called in both samples. It is then rescaled against the mean similarity of all pairs: `D = 1 - (M - MB) / (1 - MB)`.
- **Bootstrap:** 1,000 replicates. Each SNP gets a Poisson(1) weight per replicate, which approximates resampling SNPs with replacement. Genotypes are read in 20,000-SNP chunks. Setting `BLOCK_BOOTSTRAP <- TRUE` gives a block bootstrap with one weight per chunk instead.
- **Tree and support:** The NJ tree (`ape::nj`) is built from the full-data distance. Node support is the percentage of bootstrap trees that contain each clade (`ape::prop.clades`).
- **Rooting:** The tree is rooted on PI538450, IG116176, KU-1982 and TA873(Hakan). The script stops if these four don't form a single group (monophyletic) in the tree. Random seed: 48549.

## Run

```bash
sbatch run_NJ_phylogeny.sh
```

Progress is saved every 1,000,000 SNPs (`timopheevii_NJ_checkpoint.RData`). If the job is restarted, it resumes from the checkpoint and produces the same trees as an uninterrupted run. The checkpoint is deleted when the run finishes.

## Outputs

| File | Contents |
| --- | --- |
| `timopheevii_NJ_rooted_bootstrap_NJ.nwk` | Rooted NJ tree with bootstrap support (%) on nodes |
| `timopheevii_NJ_unrooted_bootstrap_NJ.nwk` | Unrooted NJ tree with bootstrap support (%) |
| `timopheevii_NJ_bootstrap_trees.nwk` | The 1,000 bootstrap trees |
| `timopheevii_NJ_distmat_ref.txt` | Full-data distance matrix |
