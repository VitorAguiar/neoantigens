#!/bin/bash

# path completo para o bam input
bam=$1

# criar uma ID para os arquivos de saida
tmpid=$(basename "$bam")
id="${tmpid%.*}"

# comandos da analise
samtools view $bam 'chr6:29222775-33629084' -b -o ${id}_mhc.bam
