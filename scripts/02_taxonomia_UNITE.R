# 02_taxonomia_UNITE.R
# Asignar una taxonomía a cada ASV a partir del resultado de BLAST (blastn de BLAST+ 2.17.0) contra la base UNITE (versión general del 19-02-2025).
# Para cada ASV se conserva el alineamiento de mayor puntaje (bitscore; si hay empate, el de menor e-value).
# La asignación se acepta solo si la identidad es de al menos 97 % y la cobertura es de al menos 60 %,
# donde la cobertura es el largo del alineamiento dividido por el largo del ASV. Si no se cumple, el ASV queda como "Unassigned".
# Entrada: datos/ASVs_vs_UNITE.tsv (salida de blastn) y datos/ASV_table.csv.
# Salida: datos/ASV_taxonomy_UNITE97_60.csv (una fila por ASV) y datos/ASV_sample_taxonomy_UNITE97_60.csv (muestra, ASV y taxonomía).

suppressPackageStartupMessages({
  library(data.table)
})

# Ruta base de la carpeta dada2_filt
filt_path <- "datos"

# Aseguramos carpeta logs
log_dir <- file.path(filt_path, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

# Comenzamos log de R
logfile <- file.path(log_dir, "R_taxonomia_UNITE.log")
sink(logfile, split = TRUE)
cat("== Iniciando anotación taxonómica con UNITE ==\n")
cat("Fecha/hora:", as.character(Sys.time()), "\n\n")

# 1) Leer tabla de "BLAST" (formato especial UNITE)
blast_file <- file.path(filt_path, "ASVs_vs_UNITE.tsv")
cat("Leyendo archivo BLAST (formato UNITE especial):", blast_file, "\n")

if (!file.exists(blast_file)) {
  sink()
  stop("No se encontró ASVs_vs_UNITE.tsv en datos (salida de blastn contra UNITE).")
}

bt <- fread(blast_file, header = FALSE)
cat("Columnas detectadas en tabla BLAST:", ncol(bt), "\n")
cat("Total de filas en tabla BLAST:", nrow(bt), "\n\n")

if (nrow(bt) == 0) {
  cat("ATENCIÓN: La tabla BLAST está vacía (0 filas). Todos los ASVs se marcarán como Unassigned.\n\n")
}

# Necesitamos al menos 9 columnas, porque toda la info útil está en V9
if (ncol(bt) < 9) {
  sink()
  stop("La tabla BLAST tiene menos de 9 columnas. El formato no coincide con el esperado de UNITE.")
}

# V1 = qseqid (ID del ASV, p.ej. "ASV0001")
# V9 = string del tipo:
# "s__Especie\tpident\tlength\tqlen\tevalue\tbitscore\tTitulo_unite"
cat("Ejemplo de V1 y V9 (primeras 3 filas):\n")
print(bt[1:3, .(V1, V9)])

# En V9 los separadores vienen como texto "\t" (barra + t).
# Primero los convertimos a tabs reales, luego partimos por tab.
v9_fixed <- gsub("\\\\t", "\t", bt[[9]])

# Separar columnas por tabulador
spl <- tstrsplit(v9_fixed, "\t")

cat("N° de columnas obtenidas al partir V9:", length(spl), "\n")

if (length(spl) < 7) {
  cat("Ejemplo de V9 original:\n")
  print(bt[1:3, V9])
  sink()
  stop("Al partir V9 por tab no se obtuvieron las 7 partes esperadas (s__, pident, length, qlen, evalue, bitscore, stitle). Revisa el formato de ASVs_vs_UNITE.tsv.")
}

species_tag <- spl[[1]]  # s__Nombre_especie (no lo usamos directo)
pident      <- as.numeric(spl[[2]])
aln_length  <- as.numeric(spl[[3]])
qlen        <- as.numeric(spl[[4]])
evalue      <- as.numeric(spl[[5]])
bitscore    <- as.numeric(spl[[6]])
stitle      <- spl[[7]]   # Título de UNITE con nombre de especie/género, etc.

cat("Resumen rápido de pident leídos desde V9:\n")
print(summary(pident))
cat("\nResumen rápido de qlen leídos desde V9:\n")
print(summary(qlen))
cat("\n")

# Construimos una tabla "blast_tab" estándar a partir de este formato
blast_tab <- data.table(
  qseqid   = bt[[1]],
  pident   = pident,
  length   = aln_length,
  qlen     = qlen,
  evalue   = evalue,
  bitscore = bitscore,
  stitle   = stitle
)

blast_tab[, coverage := length / qlen]

cat("Filas válidas (no NA) en pident:", sum(!is.na(pident)), "de", nrow(blast_tab), "\n")
cat("Filas válidas (no NA) en qlen:",   sum(!is.na(qlen)),   "de", nrow(blast_tab), "\n\n")

# 2) Quedarse con el mejor hit por ASV
cat("Seleccionando mejor hit por ASV (mayor bitscore, menor evalue)...\n")
setorder(blast_tab, qseqid, -bitscore, evalue)
blast_best <- blast_tab[!duplicated(qseqid)]
cat("ASVs únicos con al menos un hit en UNITE:", nrow(blast_best), "\n\n")

# Normalizar IDs de ASV: nos quedamos con lo que está antes del primer ";"
blast_best[, ASV := trimws(sub(";.*$", "", qseqid))]

# 3) Leer tabla de ASVs (DADA2)
asv_table_file <- file.path(filt_path, "ASV_table.csv")
cat("Leyendo tabla de ASVs:", asv_table_file, "\n")

if (!file.exists(asv_table_file)) {
  sink()
  stop("No se encontró ASV_table.csv en dada2_filt. Revisa el script de DADA2.")
}

asv_tab <- read.csv(asv_table_file, row.names = 1, check.names = FALSE)
cat("Dimensiones ASV_table (muestras x ASVs):", paste(dim(asv_tab), collapse = " x "), "\n\n")

asv_ids <- colnames(asv_tab)
cat("N° de ASVs en ASV_table:", length(asv_ids), "\n")
cat("N° de ASVs únicos en BLAST_best:", length(unique(blast_best$ASV)), "\n")
cat("N° de IDs que coinciden entre ASV_table y BLAST_best: ",
    length(intersect(asv_ids, blast_best$ASV)), "\n\n")

# 4) Construir tabla de taxonomía juntando ASVs de DADA2 con mejor hit
blast_slim <- blast_best[, .(ASV, pident, coverage, taxon_raw = stitle)]

tax_dt <- merge(
  data.table(ASV = asv_ids),
  blast_slim,
  by = "ASV",
  all.x = TRUE
)

# 5) Aplicar criterios de identidad y cobertura (paper: 97% y 60%)
cat("Aplicando criterios: pident >= 97 y coverage >= 0.60...\n")

tax_dt[, asignacion := "Unassigned"]
tax_dt[!is.na(pident) & pident >= 97 & coverage >= 0.60, asignacion := "Assigned"]

tax_dt[asignacion == "Assigned",  taxon := taxon_raw]
tax_dt[asignacion == "Unassigned", taxon := "Unassigned"]

cat("ASVs Assigned:", sum(tax_dt$asignacion == "Assigned"), "\n")
cat("ASVs Unassigned:", sum(tax_dt$asignacion == "Unassigned"), "\n\n")

# 6) Guardar tabla de taxonomía por ASV
tax_out_file <- file.path(filt_path, "ASV_taxonomy_UNITE97_60.csv")
cat("Guardando tabla de taxonomía por ASV en:\n  ", tax_out_file, "\n\n")
write.csv(tax_dt[, .(ASV, pident, coverage, taxon, asignacion)],
          tax_out_file, row.names = FALSE)

# 7) Tabla larga muestra-ASV-taxonomía
cat("Generando tabla larga muestra-ASV-taxonomía...\n")

asv_dt <- as.data.table(asv_tab, keep.rownames = "sample")
asv_long <- melt(asv_dt, id.vars = "sample",
                 variable.name = "ASV", value.name = "count")
asv_long <- asv_long[count > 0]

asv_long[tax_dt, `:=`(
  pident     = i.pident,
  coverage   = i.coverage,
  taxon      = i.taxon,
  asignacion = i.asignacion
), on = "ASV"]

asv_long_file <- file.path(filt_path, "ASV_sample_taxonomy_UNITE97_60.csv")
cat("Guardando tabla muestra-ASV-taxonomía en:\n  ", asv_long_file, "\n\n")
write.csv(asv_long, asv_long_file, row.names = FALSE)

# Respaldos en .rds
saveRDS(tax_dt,   file.path(filt_path, "ASV_taxonomy_UNITE97_60.rds"))
saveRDS(asv_long, file.path(filt_path, "ASV_sample_taxonomy_UNITE97_60.rds"))

cat("== Anotación taxonómica completada ==\n\n")
sessionInfo()
sink()
