# 01_dada2_pipeline_ITS2.R
# Obtener las variantes exactas de secuencia (ASV) a partir de las lecturas sin cebadores, con el paquete DADA2 1.30.0 (R 4.5).
# Entrada: lecturas/sample_ADTN#_R1_v2.fastq y _R2_v2.fastq (salida del paso 00).
# Salida: en datos/, filterAndTrim_counts.csv, dada2_tracking.csv, ASV_table.csv y ASVs_nochim.fasta (los usan los pasos siguientes);
# en lecturas/dada2_filt, las lecturas filtradas, los gráficos de calidad y las tablas de secuencias (.rds).
# Pasos: filtrar por calidad, estimar los errores de secuenciación, inferir las ASV, unir las lecturas pareadas,
# eliminar quimeras, contar las lecturas que quedan en cada paso y guardar la tabla de ASV y sus secuencias.

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(data.table)
})

# === Paths ===
basedir <- "lecturas"
outdir  <- "datos"
dir.create(outdir, showWarnings = FALSE)

# Archivos ya recortados solo primers/adaptadores (Cutadapt): sample_ADTN1-4_R{1,2}_v2.fastq
samples <- sprintf("sample_ADTN%d", 1:4)
fnFs <- file.path(basedir, paste0(samples, "_R1_v2.fastq"))
fnRs <- file.path(basedir, paste0(samples, "_R2_v2.fastq"))

# Chequeo de existencia
stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))

# Carpetas de salida
filt_path <- file.path(basedir, "dada2_filt")
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)

filtFs <- file.path(filt_path, paste0(samples, "_R1_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(samples, "_R2_filt.fastq.gz"))

# === (Opcional) Perfiles de calidad previos ===
# pdf(file.path(filt_path, "qual_F_raw.pdf"));  for (f in fnFs) print(plotQualityProfile(f)); dev.off()
# pdf(file.path(filt_path, "qual_R_raw.pdf"));  for (f in fnRs) print(plotQualityProfile(f)); dev.off()

# === 1) Filtrar las lecturas por calidad ===
# Se usa la función filterAndTrim. Donde:
#  - truncLen = c(0, 0): Indica que las lecturas no se cortan a un largo fijo, porque el largo de la región ITS2 cambia entre hongos.
#  - maxN = 0: Descarta las lecturas que tengan alguna base indeterminada (N).
#  - maxEE = c(2, 2): Especifica el máximo de errores esperados por lectura (2 en la forward y 2 en la reverse); las lecturas con más errores se descartan.
#  - truncQ = 11: Corta la lectura en la primera base con calidad igual o menor a 11.
#  - minLen = 50: Descarta las lecturas que quedan con menos de 50 bases después del corte.
#  - rm.phix = TRUE: Elimina las lecturas del fago PhiX, que se usa como control en la secuenciación Illumina.
#  - compress = TRUE y multithread = TRUE: Guardan los archivos comprimidos y usan varios núcleos del procesador.
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen = c(0, 0),   # sin truncación fija al inicio
                     maxN = 0,             # descarta Ns
                     maxEE = c(2, 2),      # errores esperados forward/reverse
                     truncQ = 11,          # recorte suave por calidad en extremos
                     minLen = 50,          # evita fragmentos demasiado cortos
                     rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
write.csv(out, file.path(outdir, "filterAndTrim_counts.csv"), row.names = TRUE)

# (Opcional) Perfiles de calidad tras filtro
# pdf(file.path(filt_path, "qual_F_filt.pdf"));  for (f in filtFs) print(plotQualityProfile(f)); dev.off()
# pdf(file.path(filt_path, "qual_R_filt.pdf"));  for (f in filtRs) print(plotQualityProfile(f)); dev.off()

# === 2) Estimar los errores de secuenciación e inferir las ASV ===
# learnErrors estima el modelo de errores a partir de las mismas lecturas; set.seed(1) fija la semilla para que el resultado se pueda repetir.
# derepFastq junta las lecturas idénticas y dada infiere las ASV de cada muestra.
# pool = "pseudo": Indica que en la inferencia se usa también la información de las otras muestras, lo que ayuda a detectar variantes poco abundantes.
set.seed(1)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

ddsF <- derepFastq(filtFs); names(ddsF) <- samples
ddsR <- derepFastq(filtRs); names(ddsR) <- samples

dadaF <- dada(ddsF, err = errF, multithread = TRUE, pool = "pseudo")
dadaR <- dada(ddsR, err = errR, multithread = TRUE, pool = "pseudo")

# === 3) Unir las lecturas pareadas ===
# Se usa la función mergePairs. Donde:
#  - minOverlap = 20: Exige que la lectura forward y la reverse se solapen en al menos 20 bases para unirlas.
#  - maxMismatch = 0: No permite diferencias entre las dos lecturas en la zona de solapamiento.
mergers <- mergePairs(dadaF, ddsF, dadaR, ddsR,
                      minOverlap = 20,
                      maxMismatch = 0)

# === 4) Armar la tabla de ASV y eliminar las quimeras ===
# removeBimeraDenovo elimina las quimeras (secuencias artificiales formadas en la PCR a partir de dos secuencias distintas).
# method = "consensus": Indica que una secuencia se marca como quimera según lo que ocurre en la mayoría de las muestras.
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

saveRDS(seqtab.nochim, file.path(filt_path, "seqtab_nochim.rds"))
saveRDS(seqtab,        file.path(filt_path, "seqtab_raw.rds"))

# === 5) Contar las lecturas que quedan en cada paso (dada2_tracking.csv) ===
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  input    = out[,1],
  filtered = out[,2],
  denoisedF= sapply(dadaF, getN),
  merged   = sapply(mergers, getN),
  nonchim  = rowSums(seqtab.nochim),
  row.names = samples
)
write.csv(track, file.path(outdir, "dada2_tracking.csv"))

# === 6) Ordenar las ASV de mayor a menor abundancia, nombrarlas (ASV0001, ASV0002...) y guardarlas en FASTA ===
asv_seqs <- colnames(seqtab.nochim)
asv_abun <- colSums(seqtab.nochim)
ord <- order(asv_abun, decreasing = TRUE)
asv_seqs <- asv_seqs[ord]
asv_abun <- asv_abun[ord]
asv_ids  <- sprintf("ASV%04d;size=%d;", seq_along(asv_seqs), asv_abun)

asv_fa <- DNAStringSet(asv_seqs)
names(asv_fa) <- asv_ids
Biostrings::writeXStringSet(asv_fa, filepath = file.path(outdir, "ASVs_nochim.fasta"))

# Guardar la tabla de ASV por muestra, con las columnas en el mismo orden
asv_tab <- seqtab.nochim[, ord, drop = FALSE]
colnames(asv_tab) <- sub(";size=\\d+;","", asv_ids)  # nombres cortos tipo ASV0001
write.csv(asv_tab, file.path(outdir, "ASV_table.csv"), row.names = TRUE)

message("DADA2 paso 1 completado. Tablas y secuencias en: ", outdir)
