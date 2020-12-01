# README

## Passo 1

Crie um diretório para os arquivos bam:

```console
$ mkdir bam_files
```

## Passo 2

No seu diretório de bams,

```console
$ cd bam_files
```

Crie links simbólicos para os seus Bams. Por exemplo, se os bams estiverem em 2 diretórios em lugares distintos:

```console
$ ln -s /path/to/dir1 bam_files_dir1
$ ln -s /path/to/dir2 bam_files_dir2
```

## Passo 3

Execute o script *create_bam_list.sh* para criar um arquivo-lista com todos os bams:

```console
$ ./create_bam_list.sh
```

## Passo 4

Volte para o diretório de análise:

```console
$ cd ..
```

Atualize o script *run_analysis.sh* para incluir suas análises reais.

E então execute o script com o xargs, por exemplo com 4 amostras em paralelo:

```console
$ xargs -a bam_files/bam_list.txt -n1 -P4 ./run_analysis.sh
```
