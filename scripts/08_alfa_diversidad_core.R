# =====================================================================
# 08_alfa_diversidad_core.R
# Calcular la α-diversidad de cada individuo (ADTN1 a ADTN4) y la comunidad central (core).
# Las cuatro muestras son cuatro individuos distintos de un mismo sitio (n = 4), por lo que
# los índices se resumen de forma descriptiva como media ± DE, sin pruebas estadísticas.
# Índices (paquete vegan): riqueza observada, Shannon, Simpson (1 - D), Simpson inverso y equidad de Pielou.
# También se calcula la riqueza rarefactada a la profundidad mínima, las curvas de rarefacción
# y la ocupancia: en cuántos individuos aparece cada ASV (core = presentes en los cuatro).
# Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv.
# Salidas: Alfa_diversidad_por_replica.csv y Alfa_diversidad_resumen.csv (el nombre "replica" quedó de
#   una versión anterior; cada fila es un individuo), Ocupancia_ASVs.csv (n_replicas = número de individuos),
#   Core_ASVs.csv y gráficos de trabajo en resultados/figuras/.
# =====================================================================

loc_lib <- "Rlib"
dir.create(loc_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(loc_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- c("vegan", "dplyr", "tidyr", "ggplot2", "ragg")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst)) try(install.packages(inst, lib = loc_lib, quiet = TRUE), silent = TRUE)
suppressPackageStartupMessages({ library(vegan); library(dplyr); library(tidyr); library(ggplot2) })

## 1) Rutas y estilo ---------------------------------------------
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

## 3) α-diversidad por individuo --------------------------------
reads <- rowSums(mat)
mind  <- min(reads)                       # profundidad minima (para rarefaccion)
S     <- specnumber(mat)                  # riqueza observada (nº ASVs)
H     <- diversity(mat, "shannon")        # Shannon H'
Dsi   <- diversity(mat, "simpson")        # Simpson (1 - D)
invD  <- diversity(mat, "invsimpson")     # Simpson inverso
J     <- H / log(S)                       # equidad de Pielou
Sraref<- as.numeric(rarefy(mat, mind))    # riqueza esperada a profundidad comun

alfa <- data.frame(
  Muestra = rownames(mat), Lecturas = as.integer(reads),
  Riqueza_S = as.integer(S), Riqueza_rarefactada = round(Sraref, 1),
  Shannon = round(H, 3), Simpson = round(Dsi, 4),
  InvSimpson = round(invD, 2), Pielou = round(J, 3),
  row.names = NULL)

num <- alfa[, c("Riqueza_S","Riqueza_rarefactada","Shannon","Simpson","InvSimpson","Pielou")]
resumen <- data.frame(
  Metrica = names(num),
  Media   = round(sapply(num, mean), 3),
  DE      = round(sapply(num, sd), 3),
  row.names = NULL)

write.csv(alfa,    file.path(outdir, "Alfa_diversidad_por_replica.csv"), row.names = FALSE)
write.csv(resumen, file.path(outdir, "Alfa_diversidad_resumen.csv"),     row.names = FALSE)

## 4) Comunidad central (core) y ocupancia -----------------------
occ <- colSums(mat > 0)                    # en cuántos individuos aparece cada ASV (1-4)
nrep <- nrow(mat)
nasv_by  <- sapply(1:nrep, function(k) sum(occ == k))
reads_by <- sapply(1:nrep, function(k) sum(colSums(mat)[occ == k]))
ocup <- data.frame(n_replicas = 1:nrep, n_ASVs = nasv_by, reads = reads_by,
                   pct_lecturas = round(100 * reads_by / sum(mat), 1))
write.csv(ocup, file.path(outdir, "Ocupancia_ASVs.csv"), row.names = FALSE)

core_ids <- colnames(mat)[occ == nrep]
core <- data.frame(ASV = core_ids, reads = colSums(mat)[core_ids], row.names = NULL) %>%
  left_join(tax[, c("ASV", "Genus")], by = "ASV") %>%
  arrange(desc(reads)) %>%
  mutate(pct_lecturas = round(100 * reads / sum(mat), 2))
write.csv(core, file.path(outdir, "Core_ASVs.csv"), row.names = FALSE)
core_pct <- round(100 * sum(core$reads) / sum(mat), 1)

## 5) Figuras ----------------------------------------------------
guardar <- function(p, base, w = 18, h = 11, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); print(p); dev.off()
    ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "cm", res = dpi, compression = "lzw"); print(p); dev.off()
  } else ggsave(paste0(base, ".png"), p, width = w, height = h, units = "cm", dpi = dpi)
}
guardar_base <- function(base, fun, w = 18, h = 12, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); fun(); dev.off()
    ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "cm", res = dpi, compression = "lzw"); fun(); dev.off()
  } else { png(paste0(base, ".png"), width = w, height = h, units = "cm", res = dpi); fun(); dev.off() }
}

# (a) α-diversidad: 4 índices, un punto por individuo + media±DE
long <- alfa %>% select(Muestra, Riqueza_S, Shannon, Simpson, Pielou) %>%
  pivot_longer(-Muestra, names_to = "Indice", values_to = "Valor") %>%
  mutate(Indice = factor(Indice, levels = c("Riqueza_S", "Shannon", "Simpson", "Pielou")))
p_alfa <- ggplot(long, aes(x = "A. dareoides", y = Valor)) +
  geom_jitter(width = 0.07, size = 2.2, color = "#2c7fb8") +
  stat_summary(fun = mean, fun.min = function(z) mean(z) - sd(z),
               fun.max = function(z) mean(z) + sd(z), geom = "pointrange", color = "black") +
  facet_wrap(~Indice, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "α-diversidad por individuo (media ± DE, n = 4)") +
  theme_bw(base_family = fam) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 9),
        strip.text = element_text(face = "bold"))
guardar(p_alfa, file.path(figdir, "Fig_alfa_diversidad"), w = 20, h = 9)

# (b) curva de rarefaccion (linea vertical = profundidad minima)
guardar_base(file.path(figdir, "Fig_rarefaccion"), function()
  rarecurve(mat, step = 500, sample = mind, lwd = 2,
            col = grDevices::hcl.colors(nrow(mat), "Dark3"),
            xlab = "Nº de lecturas", ylab = "ASVs (riqueza)",
            main = "Curvas de rarefacción por individuo"))

# (c) ocupancia / core
p_core <- ggplot(ocup, aes(x = factor(n_replicas), y = n_ASVs)) +
  geom_col(fill = ifelse(ocup$n_replicas == nrep, "#2c7fb8", "grey70"),
           color = "black", width = 0.72) +
  geom_text(aes(label = paste0(n_ASVs, "\n(", pct_lecturas, "%)")),
            vjust = -0.25, size = 3.1, family = fam) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Nº de individuos en que aparece el ASV",
       y = "Nº de ASVs",
       title = "Ocupancia y comunidad central (core = 4 individuos)") +
  theme_minimal(base_family = fam) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
guardar(p_core, file.path(figdir, "Fig_core_ocupancia"), w = 16, h = 11)

## 6) Reporte en consola -----------------------------------------
cat("\n================ alfa-DIVERSIDAD (n=4) ================\n")
print(alfa, row.names = FALSE)
cat("\nMedia ± DE:\n"); print(resumen, row.names = FALSE)
cat("\nProfundidad minima (rarefaccion):", mind, "lecturas\n")
cat("\n================ COMUNIDAD CENTRAL (core) ================\n")
print(ocup, row.names = FALSE)
cat("\nCORE: ", length(core_ids), " ASVs presentes en las ", nrep,
    " individuos = ", core_pct, "% de las lecturas.\n", sep = "")
cat("Generos del core (top):\n")
print(head(core %>% group_by(Genus) %>%
             summarise(ASVs = n(), reads = sum(reads), .groups = "drop") %>%
             arrange(desc(reads)), 10), row.names = FALSE)
cat("\nSalidas en:\n  ", outdir, "\n")
