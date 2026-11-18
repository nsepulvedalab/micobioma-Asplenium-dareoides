# =====================================================================
# 10_torta_asignacion.R
# Hacer la torta de lecturas asignadas y no asignadas taxonómicamente (Figura 1 de la tesis).
# Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv.
# Salida: resultados/figuras/Fig_torta_asignacion.png (y .tiff).
# =====================================================================

loc_lib <- "Rlib"
dir.create(loc_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(loc_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- c("dplyr", "tidyr", "ggplot2", "ggforce", "ragg")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst)) try(install.packages(inst, lib = loc_lib, quiet = TRUE), silent = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggforce) })

basedir <- "datos"
outdir  <- "resultados"
figdir  <- file.path(outdir, "figuras"); dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
fam <- "Calibri"
if (.Platform$OS.type == "windows") try(windowsFonts(Calibri = windowsFont("Calibri")), silent = TRUE) else fam <- ""
dec <- ","

## datos: lecturas por ASV + asignación
tax <- read.csv(file.path(basedir, "ASV_taxonomy_UNITE97_60.csv"), stringsAsFactors = FALSE, check.names = FALSE)
asv <- read.csv(file.path(basedir, "ASV_table.csv"), stringsAsFactors = FALSE, check.names = FALSE)
if (!"ASV" %in% names(tax)) names(tax)[1] <- "ASV"
names(asv)[1] <- "Sample"
reads <- asv %>% pivot_longer(-Sample, names_to = "ASV", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>% group_by(ASV) %>% summarise(reads = sum(n, na.rm = TRUE), .groups = "drop")
tax$asig <- ifelse(tax$taxon == "Unassigned" | is.na(tax$taxon) |
                   tolower(as.character(tax$asignacion)) == "unassigned",
                   "No asignados", "Asignados")
res <- reads %>% left_join(tax[, c("ASV", "asig")], by = "ASV") %>%
  mutate(asig = ifelse(is.na(asig), "No asignados", asig)) %>%
  group_by(cat = asig) %>% summarise(val = sum(reads), .groups = "drop")

## ---- geometría de etiquetas externas (igual que script 06) ----
spread_labels <- function(y, gap = 0.18, lo = -1.1, hi = 1.1) {
  n <- length(y); if (n <= 1) return(y)
  o <- order(y, decreasing = TRUE); ys <- y[o]
  for (i in 2:n) if (ys[i] > ys[i - 1] - gap) ys[i] <- ys[i - 1] - gap
  if (ys[n] < lo) { ys[n] <- lo; for (i in (n - 1):1) if (ys[i] < ys[i + 1] + gap) ys[i] <- ys[i + 1] + gap }
  out <- numeric(n); out[o] <- ys; out
}
d <- res %>% arrange(desc(val)) %>%
  mutate(cat = factor(cat, levels = cat),
         pct = 100 * val / sum(val), frac = val / sum(val),
         end = cumsum(frac) * 2 * pi, start = end - frac * 2 * pi, mid = (start + end) / 2,
         ux = sin(mid), uy = cos(mid), side = ifelse(ux >= 0, 1, -1))
d$ylab <- d$uy
for (s in c(-1, 1)) { idx <- which(d$side == s); if (length(idx)) d$ylab[idx] <- spread_labels(d$uy[idx]) }
d <- d %>% mutate(x_edge = ux, y_edge = uy, x_col = side * 1.16, x_seg2 = side * 1.30,
                  x_txt = side * 1.33, hj = ifelse(side > 0, 0, 1),
                  lab = paste0(formatC(pct, format = "f", digits = 1, decimal.mark = dec), " %"))
cols <- c("Asignados" = "#3288BD", "No asignados" = "grey80")[levels(d$cat)]

p <- ggplot(d) +
  ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, start = start, end = end, fill = cat),
                        color = "black", linewidth = 0.3) +
  geom_segment(aes(x = x_edge, y = y_edge, xend = x_col, yend = ylab), linewidth = 0.3, color = "grey25") +
  geom_segment(aes(x = x_col, y = ylab, xend = x_seg2, yend = ylab), linewidth = 0.3, color = "grey25") +
  geom_text(aes(x = x_txt, y = ylab, label = lab, hjust = hj), family = fam, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = cols, name = "ASIGNACIÓN") +
  coord_fixed(xlim = c(-1.95, 1.95), ylim = c(-1.35, 1.35), clip = "off") +
  theme_void(base_family = fam) +
  theme(legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9), plot.margin = margin(6, 10, 6, 10))

guardar <- function(p, base, w = 16, h = 12, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); print(p); dev.off()
    ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "cm", res = dpi, compression = "lzw"); print(p); dev.off()
  } else ggsave(paste0(base, ".png"), p, width = w, height = h, units = "cm", dpi = dpi)
}
guardar(p, file.path(figdir, "Fig_torta_asignacion"))
cat("\nAsignación (por lecturas):\n"); print(res %>% mutate(pct = round(100 * val / sum(val), 1)), row.names = FALSE)
cat("\nFigura guardada en:", file.path(figdir, "Fig_torta_asignacion.png"), "\n")
