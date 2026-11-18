# =====================================================================
# 09_estructura_muestras.R
# Comparar la composición de los cuatro individuos entre sí.
# - Disimilitud de Bray-Curtis (sobre abundancia relativa) y de Jaccard (presencia y ausencia).
# - Agrupamiento jerárquico por enlace promedio (dendrograma).
# - Matriz de ASV compartidos por cada par de individuos.
# - Qué par de individuos comparte los ASV que aparecen en solo dos individuos.
# - Mapa de calor de los 30 ASV más abundantes.
# Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv.
# Salidas: Disimilitud_BrayCurtis.csv, Disimilitud_Jaccard.csv, ASVs_compartidos_matriz.csv,
#          Pares_ocupancia2.csv y gráficos de trabajo en resultados/figuras/.
# =====================================================================

loc_lib <- "Rlib"
dir.create(loc_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(loc_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- c("vegan", "dplyr", "ggplot2", "ragg")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst)) try(install.packages(inst, lib = loc_lib, quiet = TRUE), silent = TRUE)
suppressPackageStartupMessages({ library(vegan); library(dplyr); library(ggplot2) })

## 1) Rutas / estilo ---------------------------------------------
basedir <- "datos"
outdir  <- "resultados"
figdir  <- file.path(outdir, "figuras")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
fam <- "Calibri"
if (.Platform$OS.type == "windows") try(windowsFonts(Calibri = windowsFont("Calibri")), silent = TRUE) else fam <- ""

## 2) Datos ------------------------------------------------------
asv <- read.csv(file.path(basedir, "ASV_table.csv"), check.names = FALSE, stringsAsFactors = FALSE)
names(asv)[1] <- "Sample"
mat <- as.matrix(asv[, -1]); rownames(mat) <- asv$Sample
storage.mode(mat) <- "double"

tax <- read.csv(file.path(basedir, "ASV_taxonomy_UNITE97_60.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!"ASV" %in% names(tax)) names(tax)[1] <- "ASV"
fp <- ifelse(tax$taxon == "Unassigned" | is.na(tax$taxon), NA_character_, sub("\\|.*", "", tax$taxon))
tax$Genus <- ifelse(is.na(fp), "Unassigned", sub("_.*", "", fp))
genus_map <- setNames(tax$Genus, tax$ASV)

## helpers de guardado -------------------------------------------
guardar <- function(p, base, w = 16, h = 18, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); print(p); dev.off()
    ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "cm", res = dpi, compression = "lzw"); print(p); dev.off()
  } else ggsave(paste0(base, ".png"), p, width = w, height = h, units = "cm", dpi = dpi)
}
guardar_base <- function(base, fun, w = 16, h = 11, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); fun(); dev.off()
    ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "cm", res = dpi, compression = "lzw"); fun(); dev.off()
  } else { png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); fun(); dev.off() }
}

## 3) Disimilitud Bray-Curtis y Jaccard --------------------------
matrel <- mat / rowSums(mat)            # abundancia relativa (corrige profundidad)
pan    <- (mat > 0) * 1                  # presencia/ausencia
bc <- as.matrix(vegdist(matrel, method = "bray"))
ja <- as.matrix(vegdist(pan,    method = "jaccard"))
write.csv(round(bc, 3), file.path(outdir, "Disimilitud_BrayCurtis.csv"))
write.csv(round(ja, 3), file.path(outdir, "Disimilitud_Jaccard.csv"))

## 4) Dendrograma ------------------------------------------------
hc <- hclust(as.dist(bc), method = "average")
guardar_base(file.path(figdir, "Fig_dendrograma_muestras"), function()
  plot(hc, hang = -1, main = "Agrupamiento de muestras (Bray-Curtis)",
       xlab = "", sub = "", ylab = "Disimilitud"))

## 5) Matriz de ASVs compartidos ---------------------------------
shared <- pan %*% t(pan)                 # diagonal = riqueza; fuera = ASVs compartidos
write.csv(shared, file.path(outdir, "ASVs_compartidos_matriz.csv"))

## 6) Tabla de pares para ASVs de ocupancia-2 --------------------
occ <- colSums(mat > 0)
asv_reads <- colSums(mat)
pair_label <- rep(NA_character_, ncol(mat))
for (j in seq_len(ncol(mat)))
  if (occ[j] == 2) pair_label[j] <- paste(sort(rownames(mat)[mat[, j] > 0]), collapse = " + ")
pair_tab <- data.frame(pair = pair_label, reads = asv_reads, stringsAsFactors = FALSE) %>%
  filter(!is.na(pair)) %>%
  group_by(pair) %>% summarise(n_ASVs = n(), reads = sum(reads), .groups = "drop") %>%
  mutate(pct_lecturas = round(100 * reads / sum(mat), 1)) %>%
  arrange(desc(reads))
write.csv(pair_tab, file.path(outdir, "Pares_ocupancia2.csv"), row.names = FALSE)

## 7) Heatmap de los ASVs dominantes -----------------------------
topN <- min(30, ncol(mat))
top_ids <- names(sort(asv_reads, decreasing = TRUE))[1:topN]
asv_ord <- hclust(dist(t(matrel[, top_ids, drop = FALSE])), "average")$order
top_ord <- top_ids[asv_ord]
lab_map <- setNames(paste0(top_ids, "  ", genus_map[top_ids]), top_ids)

hm <- expand.grid(ASV = top_ids, Sample = rownames(mat), stringsAsFactors = FALSE)
hm$rel <- mapply(function(a, s) matrel[s, a], hm$ASV, hm$Sample)
hm$Sample <- factor(hm$Sample, levels = hc$labels[hc$order])
hm$ASV    <- factor(hm$ASV, levels = top_ord, labels = lab_map[top_ord])
p_hm <- ggplot(hm, aes(Sample, ASV, fill = rel)) +
  geom_tile(color = "grey92") +
  scale_fill_viridis_c(name = "Abund.\nrelativa") +
  labs(title = paste0("Top ", topN, " ASVs por muestra (abundancia relativa)"), x = NULL, y = NULL) +
  theme_minimal(base_family = fam) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 9),
        panel.grid = element_blank())
guardar(p_hm, file.path(figdir, "Fig_heatmap_top_ASVs"), w = 17, h = 20)

## 8) Reporte e interpretación automática ------------------------
cat("\n================ ESTRUCTURA ENTRE MUESTRAS ================\n")
cat("Bray-Curtis (0 = idénticas, 1 = sin solape):\n"); print(round(bc, 3))
cat("\nOrden del dendrograma:", paste(hc$labels[hc$order], collapse = " | "), "\n")
cat("\nASVs compartidos por par (diagonal = riqueza):\n"); print(shared)
cat("\nPares que concentran los ASVs de ocupancia-2:\n"); print(pair_tab, row.names = FALSE)

top2 <- head(pair_tab, 2)
memb <- unique(unlist(strsplit(top2$pair, " \\+ ")))
if (nrow(top2) == 2 && length(memb) == nrow(mat) && sum(top2$pct_lecturas) >= 50) {
  cat("\n>> INDICIO DE ESTRUCTURA 2+2: los pares [", top2$pair[1], "] y [", top2$pair[2],
      "] son disjuntos, cubren las 4 muestras y concentran ", sum(top2$pct_lecturas),
      "% de las lecturas. Conviene revisar si corresponden a un factor (zona/fecha/lote).\n", sep = "")
} else {
  cat("\n>> SIN estructura 2+2 limpia: las lecturas de ocupancia-2 se reparten entre varios pares",
      "(variabilidad entre individuos sin agrupamiento claro).\n")
}
cat("\nSalidas en:\n  ", outdir, "\n")
