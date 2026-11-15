# ==============================================================================
# 02_bootstrap_NJ_tree.R
#
# Builds a rooted neighbour-joining (NJ) tree with bootstrap support from the
# T. timopheevii SNP set.
#
#   1. Allele-sharing similarity between every pair of samples, over SNPs
#      called in both: 1 - |g_i - g_j| / 2, rescaled against the mean
#      similarity of all pairs: D = 1 - (M - MB) / (1 - MB).
#   2. Reference NJ tree from all SNPs (each SNP counted once).
#   3. B bootstrap NJ trees. Each SNP gets a Poisson(1) weight per replicate,
#      which approximates resampling SNPs with replacement.
#   4. Node support = % of bootstrap trees containing each clade.
#   5. Tree rooted on the wild emmer (T. dicoccoides) outgroup.
#
# Genotypes are streamed in chunks, so the full genotype matrix is never held
# in memory. Progress is checkpointed and a restarted job resumes where it
# stopped, with the same random numbers as an uninterrupted run.
#
# Input  : timopheevii_WGS_SNPs.gds   (from 01_vcf_to_gds.R)
# Outputs: <prefix>_distmat_ref.txt               full-data distance matrix
#          <prefix>_bootstrap_trees.nwk           B bootstrap NJ trees
#          <prefix>_unrooted_bootstrap_NJ.nwk     NJ tree, support on nodes
#          <prefix>_rooted_bootstrap_NJ.nwk       rooted NJ tree, support on nodes
# Usage  : Rscript 02_bootstrap_NJ_tree.R   (or sbatch run_NJ_phylogeny.sh)
# ==============================================================================

suppressPackageStartupMessages({
    library(ape)
    library(SNPRelate)
    library(gdsfmt)
})

# ------------------------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------------------------
gds_file <- "timopheevii_WGS_SNPs.gds"
prefix   <- "timopheevii_NJ"

# Outgroup: the four confirmed wild emmer (T. dicoccoides) accessions
outgroup_samples <- c("PI538450", "IG116176", "KU-1982", "TA873_Hakan")

# Accessions excluded from the analysis
drop_samples <- c("IG46434", "IG116164", "IG116166")

B                   <- 1000L   # bootstrap replicates
CHUNK_SIZE          <- 20000L  # SNPs read per chunk (changing it changes the random weights)
CHECKPOINT_INTERVAL <- 50L     # save progress every 50 chunks (1,000,000 SNPs)

# FALSE: one Poisson(1) weight per SNP (used for the published tree).
# TRUE : one weight per CHUNK_SIZE block of SNPs (block bootstrap).
BLOCK_BOOTSTRAP <- FALSE

SEED <- 48549

checkpoint_file <- paste0(prefix, "_checkpoint.RData")

# ------------------------------------------------------------------------------
# 2. Samples, SNPs and sample pairs
# ------------------------------------------------------------------------------
set.seed(SEED)

if (!file.exists(gds_file)) stop("GDS file not found: ", gds_file, ". Run 01_vcf_to_gds.R first.")
gds     <- snpgdsOpen(gds_file)
all_ids <- read.gdsn(index.gdsn(gds, "sample.id"))

not_found <- setdiff(drop_samples, all_ids)
if (length(not_found) > 0)
    stop("Samples to drop not found in GDS: ", paste(not_found, collapse = ", "))

samples   <- all_ids[!(all_ids %in% drop_samples)]   # keeps the GDS order
n_samples <- length(samples)

missing_og <- setdiff(outgroup_samples, samples)
if (length(missing_og) > 0)
    stop("Outgroup samples not found in GDS: ", paste(missing_og, collapse = ", "))

# Drop SNPs that are monomorphic once the excluded samples are removed
snp_ids <- snpgdsSelectSNP(gds, sample.id = samples, autosome.only = FALSE,
                           remove.monosnp = TRUE, verbose = FALSE)
n_snps  <- length(snp_ids)
message("Samples: ", n_samples, " | SNPs: ", n_snps, " | Bootstrap replicates: ", B)

# Every sample pair (i, j) with i >= j, diagonal included
pair_idx <- which(lower.tri(matrix(0, n_samples, n_samples), diag = TRUE), arr.ind = TRUE)
n_pairs  <- nrow(pair_idx)
row_i    <- pair_idx[, 1]
col_j    <- pair_idx[, 2]

chunk_starts <- seq(1L, n_snps, by = CHUNK_SIZE)

# ------------------------------------------------------------------------------
# 3. Resume from checkpoint, or start fresh
# ------------------------------------------------------------------------------
# Per pair: sum of allele sharing (sim_*) and number of SNPs called in both (n_*).
# *_ref = reference tree (weight 1); *_boot = one column per bootstrap replicate.
first_chunk <- 1L

if (file.exists(checkpoint_file)) {
    ck <- new.env()
    load(checkpoint_file, envir = ck)
    if (!identical(ck$samples, samples) || ck$n_snps != n_snps || ck$B != B ||
        ck$CHUNK_SIZE != CHUNK_SIZE || ck$BLOCK_BOOTSTRAP != BLOCK_BOOTSTRAP)
        stop("Checkpoint does not match the current data or settings. Delete ",
             checkpoint_file, " and rerun.")
    sim_ref  <- ck$sim_ref;   n_ref  <- ck$n_ref
    sim_boot <- ck$sim_boot;  n_boot <- ck$n_boot
    assign(".Random.seed", ck$rng_state, envir = .GlobalEnv)   # continue the same random stream
    first_chunk <- ck$last_chunk + 1L
    message("Resuming from chunk ", first_chunk, " of ", length(chunk_starts))
    rm(ck)
} else {
    sim_ref  <- numeric(n_pairs)
    n_ref    <- numeric(n_pairs)
    sim_boot <- matrix(0, nrow = n_pairs, ncol = B)
    n_boot   <- matrix(0, nrow = n_pairs, ncol = B)
}

# ------------------------------------------------------------------------------
# 4. Stream genotypes and accumulate pairwise allele sharing
# ------------------------------------------------------------------------------
if (first_chunk <= length(chunk_starts)) {
    for (k in first_chunk:length(chunk_starts)) {
        s <- chunk_starts[k]
        e <- min(s + CHUNK_SIZE - 1L, n_snps)

        # samples x SNPs, genotype coded 0/1/2, NA = missing
        geno <- snpgdsGetGeno(gds, sample.id = samples, snp.id = snp_ids[s:e], verbose = FALSE)

        g_i <- geno[row_i, , drop = FALSE]
        g_j <- geno[col_j, , drop = FALSE]

        called <- (!is.na(g_i) & !is.na(g_j)) * 1          # 1 if both samples called
        diff   <- abs(g_i - g_j)
        diff[is.na(diff)] <- 0
        share  <- called * (1 - diff / 2)                   # allele sharing per pair x SNP

        share_sum  <- rowSums(share)
        called_sum <- rowSums(called)

        sim_ref <- sim_ref + share_sum
        n_ref   <- n_ref   + called_sum

        if (BLOCK_BOOTSTRAP) {
            w <- rpois(B, lambda = 1)                       # one weight per replicate for this block
            sim_boot <- sim_boot + outer(share_sum,  w)
            n_boot   <- n_boot   + outer(called_sum, w)
        } else {
            w <- matrix(rpois((e - s + 1L) * B, lambda = 1), ncol = B)   # SNPs x replicates
            sim_boot <- sim_boot + share  %*% w
            n_boot   <- n_boot   + called %*% w
        }

        rm(geno, g_i, g_j, called, diff, share, w)
        if (k %% 20 == 0) invisible(gc(verbose = FALSE))

        if (k %% CHECKPOINT_INTERVAL == 0 || e == n_snps) {
            last_chunk <- k
            rng_state  <- .Random.seed
            save(sim_ref, n_ref, sim_boot, n_boot, last_chunk, rng_state,
                 samples, n_snps, B, CHUNK_SIZE, BLOCK_BOOTSTRAP, file = checkpoint_file)
            message("Checkpoint saved at SNP ", e, " / ", n_snps)
        }
    }
}
snpgdsClose(gds)

# ------------------------------------------------------------------------------
# 5. Distance matrices and NJ trees
# ------------------------------------------------------------------------------
# Turns per-pair sums into a symmetric distance matrix
make_dist <- function(sim_sum, n_sum) {
    M  <- matrix(0, n_samples, n_samples)
    N  <- matrix(0, n_samples, n_samples)
    M[pair_idx] <- sim_sum
    N[pair_idx] <- n_sum
    M <- M + t(M);  diag(M) <- diag(M) / 2
    N <- N + t(N);  diag(N) <- diag(N) / 2
    M <- M / N
    M[N == 0] <- NA                                      # pair shares no called SNP
    MB <- mean(M, na.rm = TRUE)
    D  <- 1 - (M - MB) / (1 - MB)
    diag(D) <- 0
    D[is.na(D)] <- max(D, na.rm = TRUE) * 1.01
    dimnames(D) <- list(samples, samples)
    D
}

message("Building the reference NJ tree...")
D_ref <- make_dist(sim_ref, n_ref)
write.table(D_ref, paste0(prefix, "_distmat_ref.txt"), sep = "\t", quote = FALSE, col.names = NA)
nj_tree <- nj(as.dist(D_ref))

message("Building ", B, " bootstrap NJ trees...")
boot_trees <- lapply(seq_len(B), function(b) nj(as.dist(make_dist(sim_boot[, b], n_boot[, b]))))
class(boot_trees) <- "multiPhylo"
write.tree(boot_trees, file = paste0(prefix, "_bootstrap_trees.nwk"))

# ------------------------------------------------------------------------------
# 6. Bootstrap support
# ------------------------------------------------------------------------------
support <- prop.clades(nj_tree, boot_trees)
if (any(is.na(support))) {
    message(sum(is.na(support)), " node(s) not found in any bootstrap tree; support set to 0.")
    support[is.na(support)] <- 0
}
nj_tree$node.label <- round(support / B * 100)
write.tree(nj_tree, file = paste0(prefix, "_unrooted_bootstrap_NJ.nwk"))

# ------------------------------------------------------------------------------
# 7. Root on the outgroup
# ------------------------------------------------------------------------------
if (!is.monophyletic(nj_tree, outgroup_samples))
    stop("Outgroup is not monophyletic in the NJ tree. Inspect the unrooted tree before rooting.")

# edgelabel = TRUE keeps each support value on its branch after rerooting
rooted_tree <- root(nj_tree, outgroup = outgroup_samples, resolve.root = TRUE, edgelabel = TRUE)
rooted_file <- paste0(prefix, "_rooted_bootstrap_NJ.nwk")
write.tree(rooted_tree, file = rooted_file)

stopifnot(Ntip(read.tree(rooted_file)) == n_samples)   # file reads back with every sample

if (file.exists(checkpoint_file)) invisible(file.remove(checkpoint_file))
message("Done. Rooted tree: ", rooted_file)
