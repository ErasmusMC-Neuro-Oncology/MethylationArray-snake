#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot_CNV.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Genome-wide segmented copy number plots, one per sample, drawn from the CNV
# results stored in adata$uns by Preprocess_idat.R (section 2.8).
#
# Expects, from sesame::cnSegmentation via RUN_CNV = TRUE:
#   uns[['cnv_bin_coords_by_platform']]  list of data.frames, one per platform
#   uns[['cnv_bin_signals_by_platform']] list of data.frames, rows = samples
#   uns[['cnv_segments']]                tidy data.frame, all samples/platforms
#   uns[['cnv_failures']]                optional, sample + error
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: envs/CNV.yaml
# Usage: snakemake script directive
#
# TODO:
# 1)
#
# History:
#  18-08-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(anndata))
suppressMessages(library(ggplot2))
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input     <- snakemake@input[["adata"]]
    outfiles  <- unlist(snakemake@output[["plots"]])
    samples   <- unlist(snakemake@params[["samples"]])
    outdir    <- snakemake@params[["outdir"]]
    suffix    <- snakemake@params[["suffix"]]
    ymin      <- as.numeric(snakemake@params[["ymin"]])
    ymax      <- as.numeric(snakemake@params[["ymax"]])
    plot_w    <- as.numeric(snakemake@params[["width"]])
    plot_h    <- as.numeric(snakemake@params[["height"]])
}else{
    input     <- '/trinity/home/r115502/SSLOWGRADE/output/Methylation/methylation_data.h5ad'
    outdir    <- '/trinity/home/r115502/SSLOWGRADE/output/Methylation/CNV'
    suffix    <- '_CNV.png'
    samples   <- NULL
    outfiles  <- NULL
    ymin      <- -1.5
    ymax      <-  1.5
    plot_w    <- 14
    plot_h    <- 4
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

CHROM_ORDER <- paste0("chr", c(1:22, "X", "Y"))
#-------------------------------------------------------------------------------
# 0.3 Define functions
#-------------------------------------------------------------------------------
pick_col <- function(df, candidates, what){
    "Resolve a column by trying known aliases; sesame/DNAcopy naming varies."
    hit <- candidates[candidates %in% colnames(df)]
    if(length(hit) == 0){
        stop(sprintf("Could not find the %s column. Tried: %s. Available: %s",
                     what, paste(candidates, collapse = ", "),
                     paste(colnames(df), collapse = ", ")))
    }
    hit[1]
}

normalise_chrom <- function(x){
    x <- as.character(x)
    x <- ifelse(grepl("^chr", x), x, paste0("chr", x))
    x <- sub("^chr23$", "chrX", x)
    x <- sub("^chr24$", "chrY", x)
    x
}

placeholder_plot <- function(sample, reason){
    "Emit a labelled empty panel so the declared output always exists."
    ggplot() +
        annotate("text", x = 0, y = 0, size = 5,
                 label = sprintf("%s\n\nno copy number data\n(%s)", sample, reason)) +
        theme_void()
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
adata <- read_h5ad(input)

if(!"cnv_segments" %in% names(adata$uns)){
    stop(paste0("adata$uns has no 'cnv_segments'. Preprocess_idat.R was run ",
                "with RUN_CNV <- FALSE; rerun it with CNV segmentation on."))
}

segments     <- as.data.frame(adata$uns[["cnv_segments"]])
bin_coords   <- adata$uns[["cnv_bin_coords_by_platform"]]
bin_signals  <- adata$uns[["cnv_bin_signals_by_platform"]]
failures     <- if("cnv_failures" %in% names(adata$uns)){
                    as.data.frame(adata$uns[["cnv_failures"]])
                } else NULL

if(is.null(samples)) samples <- rownames(adata$obs)

cat(sprintf("[CNV] %d segments across %d samples; plotting %d sample(s)\n",
            nrow(segments), length(unique(segments$sample)), length(samples)))
#-------------------------------------------------------------------------------
# 1.2 Resolve column names
#-------------------------------------------------------------------------------
seg_sample <- pick_col(segments, c("sample", "ID", "Sample_Name"), "segment sample")
seg_chrom  <- pick_col(segments, c("chrom", "chrm", "seqnames", "chr",
                                   "chromosome"), "segment chromosome")
seg_start  <- pick_col(segments, c("loc.start", "start", "Start"), "segment start")
seg_end    <- pick_col(segments, c("loc.end", "end", "End"), "segment end")
seg_mean   <- pick_col(segments, c("seg.mean", "seg_mean", "seg.mean.adjusted",
                                   "mean", "signal"), "segment mean")

segments <- segments %>%
    mutate(.sample = as.character(.data[[seg_sample]]),
           .chrom  = normalise_chrom(.data[[seg_chrom]]),
           .start  = as.numeric(.data[[seg_start]]),
           .end    = as.numeric(.data[[seg_end]]),
           .mean   = as.numeric(.data[[seg_mean]])) %>%
    filter(.chrom %in% CHROM_ORDER)
#-------------------------------------------------------------------------------
# 1.3 Build a genome-wide coordinate system
#-------------------------------------------------------------------------------
# Chromosome lengths are taken from the data itself (max segment/bin end) so no
# external genome build is needed and the axis matches whatever was segmented.
coord_source <- segments %>% select(.chrom, pos = .end)

bin_tables <- list()
for(plat in names(bin_coords)){
    bc <- as.data.frame(bin_coords[[plat]])
    if(nrow(bc) == 0) next
    bc_chrom <- pick_col(bc, c("seqnames", "chrom", "chrm", "chr"), "bin chromosome")
    bc_start <- pick_col(bc, c("start", "Start"), "bin start")
    bc_end   <- pick_col(bc, c("end", "End"), "bin end")
    bin_tables[[plat]] <- data.frame(
        .chrom = normalise_chrom(bc[[bc_chrom]]),
        .start = as.numeric(bc[[bc_start]]),
        .end   = as.numeric(bc[[bc_end]]),
        stringsAsFactors = FALSE
    )
    coord_source <- bind_rows(coord_source,
                              data.frame(.chrom = bin_tables[[plat]]$.chrom,
                                         pos    = bin_tables[[plat]]$.end))
}

chrom_len <- coord_source %>%
    filter(.chrom %in% CHROM_ORDER) %>%
    group_by(.chrom) %>%
    summarise(len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    mutate(.chrom = factor(.chrom, levels = CHROM_ORDER)) %>%
    arrange(.chrom) %>%
    mutate(offset = cumsum(as.numeric(len)) - as.numeric(len),
           mid    = offset + len / 2,
           .chrom = as.character(.chrom))

genome_end <- max(chrom_len$offset + chrom_len$len)
add_offset <- function(df) left_join(df, chrom_len[, c(".chrom", "offset")],
                                     by = ".chrom")

segments <- add_offset(segments) %>%
    mutate(gstart = .start + offset, gend = .end + offset)
#-------------------------------------------------------------------------------
# 1.4 Flatten the per-platform bin signals into one lookup
#-------------------------------------------------------------------------------
# bin_signals[[platform]] has one row per sample ('sample' column) and one
# column per bin, in the same order as bin_coords[[platform]].
bins_by_sample <- list()
for(plat in names(bin_signals)){
    bs <- as.data.frame(bin_signals[[plat]])
    if(nrow(bs) == 0 || is.null(bin_tables[[plat]])) next
    samp_col <- pick_col(bs, c("sample", "rowname", "Sample_Name"), "bin sample")
    sig_cols <- setdiff(colnames(bs), samp_col)

    coords <- bin_tables[[plat]]
    if(length(sig_cols) != nrow(coords)){
        warning(sprintf(paste0("Platform %s: %d signal columns vs %d bin ",
                               "coordinates; skipping bin-level points."),
                        plat, length(sig_cols), nrow(coords)))
        next
    }
    coords <- add_offset(coords) %>%
        mutate(gmid = (.start + .end) / 2 + offset)

    for(i in seq_len(nrow(bs))){
        s <- as.character(bs[[samp_col]][i])
        bins_by_sample[[s]] <- data.frame(
            gmid   = coords$gmid,
            signal = as.numeric(unlist(bs[i, sig_cols])),
            stringsAsFactors = FALSE
        ) %>% filter(!is.na(gmid), !is.na(signal))
    }
}
#-------------------------------------------------------------------------------
# 2.1 Draw one plot per sample
#-------------------------------------------------------------------------------
for(s in samples){
    outfile <- file.path(outdir, paste0(s, suffix))
    seg_s <- segments %>% filter(.sample == paste0(s,'_tumor1'))
    if(nrow(seg_s) == 0){
        reason <- if(!is.null(failures) && s %in% failures$sample){
            "segmentation failed"
        } else {
            "sample absent from cnv_segments"
        }
        warning(sprintf("%s: %s", s, reason))
        ggsave(outfile, placeholder_plot(s, reason),
               width = plot_w, height = plot_h, dpi = 150)
        next
    }

    p <- ggplot() +
        # alternating chromosome bands
        geom_rect(data = chrom_len[seq(1, nrow(chrom_len), by = 2), ],
                  aes(xmin = offset, xmax = offset + len,
                      ymin = ymin, ymax = ymax),
                  fill = "grey95", colour = NA)

    # bin-level signal, if available for this sample
    if(!is.null(bins_by_sample[[s]]) && nrow(bins_by_sample[[s]]) > 0){
        p <- p + geom_point(data = bins_by_sample[[s]],
                            aes(x = gmid, y = pmin(pmax(signal, ymin), ymax)),
                            colour = "grey55", size = 0.15, alpha = 0.4,
                            shape = 16)
    }

    p <- p +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
        geom_vline(xintercept = c(chrom_len$offset, genome_end),
                   colour = "grey75", linewidth = 0.2) +
        # segments, coloured by direction
        geom_segment(data = seg_s,
                     aes(x = gstart, xend = gend,
                         y = pmin(pmax(.mean, ymin), ymax),
                         yend = pmin(pmax(.mean, ymin), ymax),
                         colour = ifelse(.mean > 0.1, "gain",
                                  ifelse(.mean < -0.1, "loss", "neutral"))),
                     linewidth = 1.1) +
        scale_colour_manual(values = c(gain = "#B2182B", loss = "#2166AC",
                                       neutral = "grey30"),
                            guide = "none") +
        scale_x_continuous(breaks = chrom_len$mid,
                           labels = sub("^chr", "", chrom_len$.chrom),
                           expand = c(0.005, 0)) +
        coord_cartesian(xlim = c(0, genome_end), ylim = c(ymin, ymax),
                        expand = FALSE) +
        labs(title = s, x = NULL, y = expression(log[2]~ratio)) +
        theme_classic(base_size = 11) +
        theme(axis.text.x = element_text(size = 7),
              axis.ticks.x = element_blank(),
              plot.title = element_text(face = "bold", size = 11))

    ggsave(outfile, p, width = plot_w, height = plot_h, dpi = 150)
}
#-------------------------------------------------------------------------------
# 3.1 Verify every declared output exists
#-------------------------------------------------------------------------------
# Snakemake fails the job on any missing declared output, so surface a clear
# message here instead of a MissingOutputException.
if(!is.null(outfiles)){
    absent <- outfiles[!file.exists(outfiles)]
    if(length(absent) > 0){
        stop(sprintf(paste0("%d declared output(s) not written, e.g. %s.\n",
                            "Sample names in the samplesheet probably do not ",
                            "match rownames(adata$obs): %s"),
                     length(absent), absent[1],
                     paste(head(rownames(adata$obs), 3), collapse = ", ")))
    }
}

cat(sprintf("[CNV] wrote %d plot(s) to %s\n", length(samples), outdir))
