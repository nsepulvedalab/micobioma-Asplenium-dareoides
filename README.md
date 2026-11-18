# Micobioma radicular de *Asplenium dareoides* (metabarcoding ITS2)

Datos procesados y scripts de análisis de la tesis de pregrado en Bioquímica de Nicolás A. Sepúlveda (Universidad Austral de Chile): *Caracterización del micobioma radicular y desarrollo de protocolos de establecimiento de* Asplenium dareoides*: estrategias biotecnológicas para la propagación y conservación de helechos nativos*.

## Muestras y lecturas crudas
- Raíces de cuatro individuos de *A. dareoides* (ADTN1–ADTN4), Fundo Teja Norte (Valdivia, Chile).
- Región ITS2, cebadores ITS86F (5′-GTGAATCATCGAATCTTTGAA-3′) e ITS4 (5′-TCCTCCGCTTATTGATATGC-3′); Illumina NextSeq 2000, lecturas pareadas.
- Lecturas crudas: NCBI Sequence Read Archive, BioProject **PRJNA1560947**.

## Contenido
| Carpeta | Archivo | Descripción |
|---|---|---|
| `datos/` | `ASV_table.csv` | Tabla de ASV (muestras en filas, 280 ASV en columnas; lecturas no quiméricas). |
| | `ASVs_nochim.fasta` | Secuencias representativas de los 280 ASV. |
| | `ASV_taxonomy_UNITE97_60.csv` | Asignación taxonómica por ASV (identidad, cobertura, taxón). |
| | `dada2_tracking.csv`, `filterAndTrim_counts.csv` | Retención de lecturas por etapa de DADA2. |
| `resultados/` | `*.csv` | Tablas de diversidad, ocupancia, disimilitud y anotación funcional. |
| | `ASV_FungalTraits.csv`, `Tabla_genero_FungalTraits.csv`, `Resumen_FungalTraits.csv`, `Figura3_FungalTraits.png` | Anotación funcional con FungalTraits usada en la tesis (script 13). |
| `scripts/` | ver abajo | Flujo de análisis en orden. |

## Flujo de análisis
| Paso | Script | Qué hace |
|---|---|---|
| 0 | `00_cutadapt.bat` | Remueve los cebadores con Cutadapt 5.1 (Python 3.11.5): solapamiento mínimo 10 pb, error máximo 0,1. |
| 1 | `01_dada2_pipeline_ITS2.R` | DADA2 1.30.0 (R 4.5): filtrado sin truncamiento fijo (maxEE = 2, truncQ = 11, minLen = 50, sin N, sin PhiX), inferencia con `pool = "pseudo"`, fusión con solapamiento ≥ 20 pb sin diferencias y remoción de quimeras por consenso. Guarda la tabla de ASV, las secuencias y el seguimiento de lecturas en `datos/`. |
| 1b | `01b_blast_UNITE.bat` | `blastn` (BLAST+ 2.17.0) de `ASVs_nochim.fasta` contra UNITE (versión general del 19-02-2025); salida tabular `-outfmt "6 qseqid sseqid pident length qlen evalue bitscore stitle"` en `datos/ASVs_vs_UNITE.tsv`. |
| 2 | `02_taxonomia_UNITE.R` | Mejor alineamiento por ASV (bitscore); umbrales de identidad ≥ 97 % y cobertura ≥ 60 % (longitud del alineamiento / longitud del ASV). |
| 6 | `06_figuras_tortas.R` | Tortas de composición por género y especie. |
| 8 | `08_alfa_diversidad_core.R` | Riqueza, Shannon, Simpson, Pielou (vegan 2.7-5), rarefacción y comunidad central. |
| 9 | `09_estructura_muestras.R` | Bray-Curtis, Jaccard, agrupamiento jerárquico y mapa de calor. |
| 10 | `10_torta_asignacion.R` | Proporción de lecturas asignadas y no asignadas. |
| 11 | `11_figura2_composicion.py` | Figura 2 de la tesis (Python 3, pandas, matplotlib). |
| 12 | `12_figuras_4a8.py` | Figuras 4 a 8 de la tesis: α-diversidad, rarefacción, ocupancia, dendrograma de Bray-Curtis y mapa de calor (Python 3, numpy, scipy, matplotlib). |
| 13 | `13_funcional_FungalTraits.py` | Anotación funcional por género con FungalTraits (Põlme et al., 2020; Tabla S1 «Traits of genera»): estilo de vida primario y secundario, capacidad endofítica y capacidad patogénica para plantas; modo trófico agrupado desde el estilo de vida primario (saprótrofo, patótrofo, simbiótrofo); proporciones ponderadas por lecturas. Genera la Figura 3 de la tesis (Python 3, pandas, openpyxl, matplotlib). |

Los pasos 3, 4, 5 y 7 del flujo original no se incluyen porque no forman parte de los resultados de la tesis.

Los scripts en R (06 y 08–10) calculan las tablas de `resultados/` y generan gráficos de trabajo. Las figuras de la tesis se generan con los scripts en Python 11 y 12, que repiten los mismos cálculos (las cifras coinciden con las tablas de R) con formato en castellano. La Figura 1 de la tesis proviene del script 10 y la Figura 3, del script 13.

## Uso
1. Descargar las lecturas del BioProject PRJNA1560947 a una carpeta `lecturas/` junto a `scripts/` (nombres `sample_ADTN#_R1.fastq` / `_R2.fastq`, descomprimidos).
2. Ejecutar los scripts en orden desde la carpeta raíz del repositorio. Las rutas son relativas (`datos/`, `resultados/`, `lecturas/`, `Rlib/`).
3. UNITE y FungalTraits no se redistribuyen: se descargan desde sus fuentes oficiales. Para el paso 1b, descargar UNITE (versión general del 19-02-2025, https://unite.ut.ee) y dejar el archivo `UNITE_sh_general_19.02.2025.fasta` en `datos_externos/`. Para el script 13, descargar la Tabla S1 «Traits of genera» del material suplementario de Põlme et al. (2020) (https://doi.org/10.1007/s13225-020-00466-2) y dejarla como `datos_externos/13225_2020_466_MOESM4_ESM.xlsx`.

## Reproducibilidad
Al repetir el flujo completo desde las lecturas crudas (pasos 0 a 2) se obtuvieron los mismos 280 ASV, con las mismas secuencias y la misma asignación taxonómica. Los conteos por ASV pueden variar en unas pocas lecturas entre ejecuciones de DADA2 (en la prueba, 538.036 frente a 538.025 lecturas en total); las tablas de `datos/` corresponden a la ejecución usada en la tesis.

## Referencias
- Callahan, B.J. et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.
- Nilsson, R.H. et al. (2019). The UNITE database for molecular identification of fungi. *Nucleic Acids Research*, 47, D259–D264.
- Garnica, S., Soto-Rauch, G., Leffler, E.M., Núñez, C., Gómez-Espinoza, J., Romero, E., Ahumada-Rudolph, R. y Cabrera-Pardo, J.R. (2025). High diversity of fungal ecological groups from ice-free pristine and disturbed areas in the Fildes Peninsula, King George Island, Antarctica. PLoS One, 20, e0317571. (flujo de procesamiento de referencia).
- Põlme, S. et al. (2020). FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. *Fungal Diversity*, 105, 1–16.
