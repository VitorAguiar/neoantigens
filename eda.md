Alelos HLA em câncer de bexiga
================

``` r
library(tidyverse)
library(rvest)
```

## Alelos dos pacientes

``` r
# ler com read.table ao invés de read_tsv
# porque o arquivo tem inconsistências com o separador nas últimas linhas
# espaços em número variável ao invés de tabs
hlatypes <- read.table("./hlatypes.txt", header = TRUE) %>%
    as_tibble() %>%
    pivot_longer(-1, names_to = "locus", values_to = "allele") %>%
    mutate(locus = sub("^(HLA)\\.([ABC])_\\d$", "\\1-\\2", locus),
           allele = toupper(allele),
           allele = sub("^HLA_([^_]+)_(\\d+)_(\\d+)$", "\\1*\\2:\\3", allele),
           allele = recode(allele, "A*07:07" = "A*02:07")) #fix typo

hlatypes
```

    # A tibble: 234 x 3
       Patient            locus allele 
       <chr>              <chr> <chr>  
     1 SRR650130_B100.sra HLA-A A*31:01
     2 SRR650130_B100.sra HLA-A A*24:02
     3 SRR650130_B100.sra HLA-B B*07:02
     4 SRR650130_B100.sra HLA-B B*15:25
     5 SRR650130_B100.sra HLA-C C*07:02
     6 SRR650130_B100.sra HLA-C C*03:04
     7 SRR650131_B101.sra HLA-A A*02:03
     8 SRR650131_B101.sra HLA-A A*11:01
     9 SRR650131_B101.sra HLA-B B*58:01
    10 SRR650131_B101.sra HLA-B B*40:01
    # … with 224 more rows

## Busca dos alelos nos bancos do NetMHC e NetMHCpan

Vamos comparar os alelos dos pacientes com aqueles presentes nos bancos
do NetMHC e NetMHCpan:

``` r
netmhc <- 
    "https://services.healthtech.dtu.dk/services/NetMHC-4.0/MHC_allele_names.txt" %>%
    read_tsv(col_names = FALSE) %>%
    setNames(c("id", "allele", "locus")) %>%
    select(locus, allele) %>%
    mutate(allele = sub("^HLA-", "", allele))

netmhcpan <- 
    "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/MHC_allele_names.txt" %>%
    read_tsv(col_names = FALSE) %>%
    setNames("allele") %>%
    mutate(allele = sub("^HLA-([A-Z]+)(.+)$", "\\1*\\2", allele))
   
hlatypes %>%
    drop_na() %>%
    mutate(NetMHC = allele %in% netmhc$allele,
           NetMHCpan = allele %in% netmhcpan$allele) %>%
    summarise_at(vars(4:5), mean) %>%
    pivot_longer(1:2, names_to = "database", values_to = "percent") %>%
    mutate(percent = scales::percent(percent)) %>%
    knitr::kable()
```

| database  | percent |
| :-------- | :------ |
| NetMHC    | 65%     |
| NetMHCpan | 100%    |

Vemos que 65% dos alelos são encontrados no NetMHC, mas todos os alelos
são encontrados no NetMHCpan v4.1.

## Busca no allelefrequencies.net para obter a frequência alélica

Vamos buscar esses alelos no Allelefrequencies.net para ver sua
frequência na população chinesa.

Para isso criei a função abaixo:

``` r
get_frequency <- function(allele, country = "") {
    
    hlaurl <- 
        "http://www.allelefrequencies.net/hla6006a.asp?hla_selection=%s&hla_country=%s"

    hlahtml <- hlaurl %>%
        sprintf(sub("HLA-", "", allele), country) %>%
        read_html()

    hlatbl <- html_node(hlahtml, "table.tblNormal")
    
    if (class(hlatbl) == "xml_missing") return(tibble(allele = allele, f = NA))
    
    hlatbl %>%
        html_table(fill = TRUE, header = TRUE) %>%
        select(2, 4, 6, 8) %>%
        as_tibble() %>%
        setNames(c("allele", "population", "f", "n")) %>%
        mutate(n = as.integer(sub(",", "", n)))
}
```

Ela tem o seguinte output quando aplicada a um alelo:

``` r
get_frequency("A*01:01", "Brazil")
```

    # A tibble: 4 x 4
      allele  population                           f     n
      <chr>   <chr>                            <dbl> <int>
    1 A*01:01 Brazil  Puyanawa                 0.043   150
    2 A*01:01 Brazil Belo Horizonte Caucasian  0.079    95
    3 A*01:01 Brazil Mixed                     0.091   108
    4 A*01:01 Brazil Vale do Ribeira Quilombos 0       144

Vamos aplicar a todos os alelos dos pacientes, para obter suas
frequências na China. Como são várias amostragens no país, vou tirar
uma média ponderada da frequência alélica dados os tamanhos amostrais,
para cada alelo.

Aplicar a todos os alelos pode levar alguns minutos, então vou
paralelizar com `furrr`.

``` r
library(furrr)

plan(multisession, workers = 8)

unique_df <- hlatypes %>% 
    drop_na(allele) %>%
    distinct(locus, allele) %>%
    arrange(locus, allele)

allele_freqs <- unique_df$allele %>%
    future_map_dfr(get_frequency, country = "China") %>%
    left_join(unique_df, .) %>%
    group_by(locus, allele) %>%
    summarise(wf = weighted.mean(f, n)) %>%
    ungroup()

#saveRDS(allele_freqs, "allele_freqs.rds")
```

``` r
allele_freqs
```

    # A tibble: 64 x 3
       locus allele         wf
       <chr> <chr>       <dbl>
     1 HLA-A A*01:01   0.0284 
     2 HLA-A A*02:01   0.112  
     3 HLA-A A*02:03   0.0367 
     4 HLA-A A*02:06   0.0455 
     5 HLA-A A*02:07   0.0936 
     6 HLA-A A*02:133 NA      
     7 HLA-A A*03:01   0.0243 
     8 HLA-A A*11:01   0.227  
     9 HLA-A A*23:01   0.00402
    10 HLA-A A*24:02   0.152  
    # … with 54 more rows

Verificamos quais alelos não estão no banco do allelefrequencies.net:

``` r
allele_freqs %>% 
    filter(is.na(wf)) %>%
    pull(allele)
```

    [1] "A*02:133"

Podemos verificar se alelo ocorre em alguma outra população (não
especificar nenhum valor pro argumento `country`):

``` r
get_frequency("A*02:133") %>%
    filter(f > 0) %>%
    mutate(f = format(f, scientific = FALSE))
```

    # A tibble: 2 x 4
      allele   population                      f           n
      <chr>    <chr>                           <chr>   <int>
    1 A*02:133 Germany DKMS - Austria minority 0.00030  1698
    2 A*02:133 Germany pop 8                   0.00001 39689

Vemos que o `A*02:133` não foi descrito na China (nesse banco de dados),
apenas possui uma frequência muita baixa numa outra população.

No entanto, para os outros alelos, podemos plotar a frequência:

``` r
library(ggthemes)

allele_freqs %>%
    mutate(allele = fct_reorder(allele, -wf)) %>%
    ggplot(aes(x = wf, y = allele, fill = wf)) +
    geom_col(color = "grey50", size = .25) +
    scale_y_discrete(limits = rev) +
    facet_wrap(~locus, scales = "free", nrow = 1) +
    scale_fill_continuous_tableau(guide = guide_colourbar(direction = "horizontal",
                                                          barwidth = 10,
                                                          barheight = .5)) +
    labs(x = "Allele frequency in China", y = NULL, fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
```

![](eda_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
