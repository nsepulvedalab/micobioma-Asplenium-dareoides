@echo off
REM 01b_blast_UNITE.bat
REM Comparar las secuencias de los ASV contra la base UNITE con BLAST+ (blastn), antes de asignar la taxonomia (paso 02).
REM Se usa la version general de UNITE del 19-02-2025 (archivo UNITE_sh_general_19.02.2025.fasta), que se descarga desde
REM https://unite.ut.ee y se deja en la carpeta datos_externos. Los comandos son:
REM   makeblastdb -in (UNITE.fasta) -dbtype nucl -out (base)
REM   blastn -query datos\ASVs_nochim.fasta -db (base) -outfmt "6 qseqid sseqid pident length qlen evalue bitscore stitle" -out datos\ASVs_vs_UNITE.tsv
REM Donde:
REM  - makeblastdb: Prepara la base UNITE en el formato que usa BLAST (se hace una sola vez).
REM  - -dbtype nucl: Indica que la base es de nucleotidos.
REM  - blastn: Compara secuencias de nucleotidos contra la base.
REM  - -query: Archivo con las secuencias de los 280 ASV (salida del paso 01).
REM  - -outfmt "6 ...": Salida en tabla, una fila por alineamiento, con estas columnas: ASV, secuencia de UNITE,
REM    porcentaje de identidad, largo del alineamiento, largo del ASV, e-value, puntaje (bitscore) y nombre completo
REM    de la secuencia de UNITE (con su taxonomia).
REM  - Se usan las demas opciones por defecto de blastn (hasta 500 secuencias de UNITE por ASV).
REM El paso 02 toma de esta tabla el mejor alineamiento de cada ASV y aplica los umbrales de 97 % de identidad y 60 % de cobertura.

cd /d "%~dp0.."

if not exist datos_externos\UNITE_db.nsq (
    makeblastdb -in datos_externos\UNITE_sh_general_19.02.2025.fasta -dbtype nucl -out datos_externos\UNITE_db
)

blastn -query datos\ASVs_nochim.fasta -db datos_externos\UNITE_db ^
    -outfmt "6 qseqid sseqid pident length qlen evalue bitscore stitle" ^
    -out datos\ASVs_vs_UNITE.tsv

echo Listo: datos\ASVs_vs_UNITE.tsv
pause
