# =====================================================================
# 06_figuras_tortas.R
# Gráficos de torta de trabajo con la composición por género y por especie.
# Estos gráficos sirvieron para revisar los resultados; la Figura 2 de la tesis
# se hizo con el script 11, con los mismos números.
# Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv.
# Salida: resultados/figuras/ (PNG y TIFF).
# =====================================================================

loc_lib <- "Rlib"
dir.create(loc_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(loc_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

need <- c("dplyr", "tidyr", "ggplot2", "ggforce", "ragg", "patchwork")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst)) try(install.packages(inst, lib = loc_lib, quiet = TRUE), silent = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggforce)
})

## ---------------------------------------------------------------
## 1) Rutas y fuente  (EDITAR ESTAS LÍNEAS SI CAMBIA LA UBICACIÓN)
## ---------------------------------------------------------------
# Carpeta con los datos (de aquí se LEE):
basedir <- "datos"
# Carpeta de salida en la carpeta conectada (aquí se ESCRIBE; revisable):
outdir  <- "resultados"
figdir  <- file.path(outdir, "figuras")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

fam <- "Calibri"
if (.Platform$OS.type == "windows") {
  try(windowsFonts(Calibri = windowsFont("Calibri")), silent = TRUE)
} else {
  fam <- ""   # Linux/Mac: usa la fuente por defecto (Calibri no existe)
}

dec <- ","    # separador decimal de los porcentajes ("," estilo español)

## ---------------------------------------------------------------
## 2) Lectura de datos base
## ---------------------------------------------------------------
tax <- read.csv(file.path(basedir, "ASV_taxonomy_UNITE97_60.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
asv <- read.csv(file.path(basedir, "ASV_table.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
if (!"ASV" %in% names(tax)) names(tax)[1] <- "ASV"
names(asv)[1] <- "Sample"

# género / especie desde 'taxon'
first_part  <- ifelse(tax$taxon == "Unassigned" | is.na(tax$taxon),
                      NA_character_, sub("\\|.*", "", tax$taxon))
tax$Genus   <- ifelse(is.na(first_part), "Unassigned", sub("_.*", "", first_part))
tax$Especie <- ifelse(is.na(first_part), "Unassigned", first_part)

# abundancia total por ASV (suma de los 4 individuos)
per_asv <- asv %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  group_by(ASV) %>% summarise(reads = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(tax[, c("ASV", "Genus", "Especie")], by = "ASV")

res_genero  <- per_asv %>% group_by(cat = Genus)   %>% summarise(val = sum(reads), .groups = "drop")
res_especie <- per_asv %>% group_by(cat = Especie) %>% summarise(val = sum(reads), .groups = "drop")

## ---------------------------------------------------------------
## 3) FUNCIONES DE TORTA
## ---------------------------------------------------------------

# Separa verticalmente las etiquetas de un mismo lado para que NO se
# solapen, conservando su orden (la de arriba sigue arriba). Empuja hacia
# abajo y, si se desborda, reacomoda hacia arriba.
spread_labels <- function(y, gap = 0.135, lo = -1.22, hi = 1.22) {
  n <- length(y); if (n <= 1) return(y)
  ord <- order(y, decreasing = TRUE)
  ys  <- y[ord]
  for (i in 2:n) if (ys[i] > ys[i - 1] - gap) ys[i] <- ys[i - 1] - gap
  if (ys[n] < lo) {
    ys[n] <- lo
    for (i in (n - 1):1) if (ys[i] < ys[i + 1] + gap) ys[i] <- ys[i + 1] + gap
  }
  if (ys[1] > hi) ys[1] <- hi
  out <- numeric(n); out[ord] <- ys; out
}

# Prepara: agrupa a top_n + "Otros", calcula % y geometría de etiquetas
# externas (ángulo medio de cada rebanada + posición vertical separada).
prep_pie <- function(df, top_n = 10, drop = NULL, gap = 0.135) {
  d <- df %>% group_by(cat) %>% summarise(val = sum(val), .groups = "drop")
  if (!is.null(drop)) d <- d %>% filter(!cat %in% drop)
  d <- d %>% arrange(desc(val))
  if (nrow(d) > top_n) {
    d <- bind_rows(d[seq_len(top_n), ],
                   tibble(cat = "Otros", val = sum(d$val[-seq_len(top_n)])))
  }
  d <- d %>%
    mutate(cat_disp = gsub("_", " ", cat),
           cat_disp = factor(cat_disp, levels = gsub("_", " ", cat)),  # conserva el orden
           pct  = 100 * val / sum(val),
           frac = val / sum(val),
           end   = cumsum(frac) * 2 * pi,
           start = end - frac * 2 * pi,
           mid   = (start + end) / 2,
           ux = sin(mid), uy = cos(mid),
           side = ifelse(ux >= 0, 1, -1))
  # separación vertical de etiquetas, por lado (izquierdo / derecho)
  d$ylab <- d$uy
  for (s in c(-1, 1)) {
    idx <- which(d$side == s)
    if (length(idx)) d$ylab[idx] <- spread_labels(d$uy[idx], gap = gap)
  }
  d <- d %>%
    mutate(x_edge = ux * 1.00, y_edge = uy * 1.00,  # punto sobre el borde de la torta
           x_col  = side * 1.16,                     # columna del codo
           x_seg2 = side * 1.30,                     # fin de la guía horizontal
           x_txt  = side * 1.33,                     # ancla del texto
           hj  = ifelse(side > 0, 0, 1),
           lab = paste0(formatC(pct, format = "f", digits = 1, decimal.mark = dec), " %"))
  d
}

# Paleta: color a cada categoría; grises DISTINTOS a Unassigned y Otros.
pal_fun <- function(levs) {
  base_levs <- setdiff(levs, c("Unassigned", "Otros", "No asignado"))
  cols <- character(0)
  if (length(base_levs))
    cols <- setNames(grDevices::hcl.colors(length(base_levs), "Spectral"), base_levs)
  if ("Unassigned"  %in% levs) cols["Unassigned"]  <- "grey80"
  if ("Otros"       %in% levs) cols["Otros"]       <- "grey55"
  if ("No asignado" %in% levs) cols["No asignado"] <- "grey80"
  cols[levs]
}

# Dibuja la torta con etiquetas externas + líneas guía en codo.
plot_pie <- function(d, legend_title, title = "", min_pct_label = 0) {
  dl <- d %>% filter(pct >= min_pct_label)   # etiquetas a mostrar
  cols <- pal_fun(levels(d$cat_disp))
  ggplot(d) +
    ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                              start = start, end = end, fill = cat_disp),
                          color = "black", linewidth = 0.3) +
    # línea guía en codo: borde de la rebanada -> columna -> etiqueta
    geom_segment(data = dl, aes(x = x_edge, y = y_edge, xend = x_col, yend = ylab),
                 linewidth = 0.3, color = "grey25") +
    geom_segment(data = dl, aes(x = x_col, y = ylab, xend = x_seg2, yend = ylab),
                 linewidth = 0.3, color = "grey25") +
    geom_text(data = dl, aes(x = x_txt, y = ylab, label = lab, hjust = hj),
              family = fam, fontface = "bold", size = 3.4) +
    scale_fill_manual(values = cols, name = toupper(legend_title)) +
    coord_fixed(xlim = c(-1.95, 1.95), ylim = c(-1.35, 1.35), clip = "off") +
    labs(title = title) +
    theme_void(base_family = fam) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
          legend.title = element_text(face = "bold", size = 10),
          legend.text  = element_text(size = 9),
          legend.key.size = unit(0.45, "cm"),
          plot.margin  = margin(6, 10, 6, 10))
}

# Exporta PNG + TIFF (usa ragg si está disponible).
guardar <- function(p, base, w = 19, h = 13, dpi = 600) {
  png_f <- paste0(base, ".png"); tif_f <- paste0(base, ".tiff")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_f, width = w, height = h, units = "cm", res = dpi)
    print(p); dev.off()
    ragg::agg_tiff(tif_f, width = w, height = h, units = "cm", res = dpi,
                   compression = "lzw"); print(p); dev.off()
  } else {
    ggsave(png_f, p, width = w, height = h, units = "cm", dpi = dpi)
    grDevices::tiff(tif_f, width = w, height = h, units = "cm", res = dpi,
                    compression = "lzw",
                    type = if (.Platform$OS.type == "windows") "cairo" else "cairo")
    print(p); dev.off()
  }
  message("Guardado: ", png_f, "  y  ", tif_f)
}

## ---------------------------------------------------------------
## 4) TORTA 1 y 2: género y especie  (incluyen 'Unassigned')
## ---------------------------------------------------------------
d_gen <- prep_pie(res_genero,  top_n = 10)
d_esp <- prep_pie(res_especie, top_n = 10)

p_gen <- plot_pie(d_gen, legend_title = "Género",  title = "Composición por género")
p_esp <- plot_pie(d_esp, legend_title = "Especie", title = "Composición por especie")

guardar(p_gen, file.path(figdir, "Fig_torta_genero"))
guardar(p_esp, file.path(figdir, "Fig_torta_especie"))

# Panel combinado A/B (si patchwork está instalado)
if (requireNamespace("patchwork", quietly = TRUE)) {
  comb <- patchwork::wrap_plots(p_gen + labs(title = NULL),
                                p_esp + labs(title = NULL), nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A")
  guardar(comb, file.path(figdir, "Fig_torta_genero_especie"), w = 34, h = 14)
}

cat("\nFiguras guardadas en:\n  ", figdir, "\n")

## ---------------------------------------------------------------
## Leyenda sugerida para la figura combinada (tesis):
## "Secuencias ITS fúngicas de raíces de Asplenium dareoides. ASVs
##  inferidos con DADA2 y anotados contra UNITE mediante BLASTn
##  (>=97 % identidad, >=60 % cobertura). (A) Composición porcentual
##  por género; (B) composición porcentual por especie."
## ---------------------------------------------------------------
