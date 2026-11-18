@echo off
REM 00_cutadapt.bat
REM Eliminar los cebadores del inicio de las lecturas de las cuatro muestras (ADTN1 a ADTN4).
REM Se usa el programa Cutadapt 5.1 (Python 3.11.5) con el siguiente comando, repetido para cada muestra:
REM   python -m cutadapt -g ^GTGAATCATCGAATCTTTGAA -G ^TCCTCCGCTTATTGATATGC -O 10 -e 0.1 --pair-filter=any -o (salida R1) -p (salida R2) (entrada R1) (entrada R2)
REM Donde:
REM  - python -m cutadapt: Ejecuta Cutadapt desde Python.
REM  - -g ^GTGAATCATCGAATCTTTGAA: Especifica el cebador forward ITS86F que se elimina de la lectura R1. El simbolo ^ indica que el cebador debe estar al inicio (extremo 5') de la lectura.
REM  - -G ^TCCTCCGCTTATTGATATGC: Especifica el cebador reverse ITS4 que se elimina de la lectura R2, tambien al inicio.
REM  - -O 10: Define el solapamiento minimo de 10 bases entre la lectura y el cebador para considerar que el cebador esta presente.
REM  - -e 0.1: Define la tasa de error maxima permitida al reconocer el cebador (10 %%).
REM  - --pair-filter=any: Indica que, si se aplica un filtro de descarte, basta con que una de las dos lecturas no lo cumpla para descartar el par.
REM  - -o y -p: Especifican los archivos de salida para R1 y R2 (sample_ADTN#_R1_v2.fastq y sample_ADTN#_R2_v2.fastq).
REM  - Los dos ultimos archivos son las lecturas crudas de cada muestra (sample_ADTN#_R1.fastq y sample_ADTN#_R2.fastq).
REM El registro de cada muestra queda en la carpeta logs.

cd /d "%~dp0..\lecturas"

REM === Crear carpeta de logs si no existe ===
if not exist logs mkdir logs

REM === Loop para procesar cada muestra ===
for %%A in (1 2 3 4) do (
    echo Procesando muestra ADTN%%A ...

    python -m cutadapt ^
        -g ^GTGAATCATCGAATCTTTGAA ^
        -G ^TCCTCCGCTTATTGATATGC ^
        -O 10 -e 0.1 --pair-filter=any ^
        -o sample_ADTN%%A_R1_v2.fastq -p sample_ADTN%%A_R2_v2.fastq ^
        sample_ADTN%%A_R1.fastq sample_ADTN%%A_R2.fastq ^
        > "logs\cutadapt_ADTN%%A.log" 2>&1
)

echo ===========================================================
echo Proceso terminado. Los archivos *_v2.fastq quedaron en la carpeta lecturas.
echo Los logs están en la carpeta "logs".
echo ===========================================================
pause
