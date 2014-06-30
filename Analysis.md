






Blau syndrome SNP sets analysis
========================================================

Analysis of individual sets
----------------------------
We analyze 3 sets of SNPs associated with NOD2 pathway in EA and AA populations:
1) `TAB1-AA` - SNPs associated with TAB1 region in AA population.
2) `TAB2-EA` - SNPs associated with TAB2 region in EA population
3) `TAB12-combined` - SNPs associated with Blau syndrome in both populations

The analysis is organized as follows:
1) Each set is tested for enrichmed co-localization in 4498 genome annotation (epigenomic) marks from the Encode project. Top 10 most significant enriched/depleted associations are shown. Look for `1Encode` postfix in the headers.
2) Each set is tested vs. 267 other genome annotation marks. Look for `2other` postfix
3) Each set is further tested for enrichment in custom categories of genomic features. Look for the description of their postfixes in the [FAQ](https://mdozmorov.github.io/grdocs/misc/Faq.html)

P-values are corrected for multiple testing using FDR. For each analysis, only the significant enrichment results (p < 0.01) are shown. Negative p-value indicates depletion. Tables show top 10 or less enriched/depleted associations, barplots visualize the same results. If an analysis does not have any significant results, it is not displayed at all.

TAB12-combined analysis
===========================
This is the most interesting analysis because it contrasts the two conditions. Look at the differences between AA and EA populations.


```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-1Encode analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 673"
##                                             Row.names pval.TAB1.AA pval.TAB2.EA                                                                    V2
## 1             wgEncodeBroadHistoneH1hescH3k36me3StdPk     4.52e-35    -1.16e-60     H1-hESC H3K36me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2                wgEncodeBroadHistoneHelas3Pol2bStdPk     2.48e-33   -1.95e-226         HeLa-S3 Pol2 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3               wgEncodeBroadHistoneHsmmH3k27me3StdPk    1.85e-145     1.02e-42        HSMM H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4              wgEncodeBroadHistoneHuvecH3k27me3StdPk    2.78e-106    6.04e-213       HUVEC H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5                 wgEncodeBroadHistoneOsteoH3k27me3Pk    4.43e-104    3.75e-204 Osteoblasts H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6              wgEncodeHaibGenotypeMyometrRegionsRep1     8.15e-29    3.18e-111             Myometr Copy number variants Replicate 1 from ENCODE/HAIB
## 7    wgEncodeSunyAlbanyGeneStGm12878Pabpc1RbpAssocRna     3.71e-76    -5.40e-89                                                                  <NA>
## 8  wgEncodeSunyAlbanyGeneStGm12878Pabpc1RbpAssocRnaV2     3.71e-76    -5.40e-89 GM12878 PABPC1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## 9      wgEncodeSunyAlbanyGeneStHelas3T7tagRbpAssocRna     3.47e-85     7.30e-27                                                                  <NA>
## 10   wgEncodeSunyAlbanyGeneStHelas3T7tagRbpAssocRnaV2     3.47e-85     7.30e-27  HeLa-S3 T7Tag RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB121.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 534"
##                                              Row.names pval.TAB1.AA pval.TAB2.EA                                                                       V2
## 1               wgEncodeBroadHistoneK562Ezh239875StdPk    -9.47e-57    6.35e-113       K562 EZH2 (39875) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2                       wgEncodeBroadHistoneK562NcorPk    -4.50e-53    3.07e-146               K562 NCoR Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3                    wgEncodeBroadHistoneK562P300StdPk    -8.70e-59    7.65e-113               K562 P300 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4                       wgEncodeBroadHistoneK562PcafPk    -8.70e-59     7.17e-99               K562 PCAF Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5           wgEncodeGisChiaPetMcf7Pol2InteractionsRep2    -2.42e-78   -2.25e-303              MCF-7 Pol2 ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 6    wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRna   -3.87e-102   -4.65e-173                                                                     <NA>
## 7  wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRnaV2   -3.87e-102   -4.65e-173 HeLa-S3 RIP-Input RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## 8        wgEncodeSunyAlbanyGeneStHepg2T7tagRbpAssocRna    6.31e-104    -1.16e-50                                                                     <NA>
## 9        wgEncodeSunyAlbanyGeneStK562Pabpc1RbpAssocRna     8.57e-91    -1.86e-61                                                                     <NA>
## 10     wgEncodeSunyAlbanyGeneStK562Pabpc1RbpAssocRnaV2     8.57e-91    -1.86e-61       K562 PABPC1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB122.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-2other analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 67"
##        Row.names pval.TAB1.AA pval.TAB2.EA                                                                                    V2
## 1        affyU95     8.67e-16    -1.87e-99                              Alignments of Affymetrix Consensus/Exemplars from HG-U95
## 2       ccdsGene     1.12e-24   -2.44e-113                                                                         Consensus CDS
## 3  coriellDelDup     8.23e-12     1.70e-33                                                Coriell Cell Line Copy Number Variants
## 4      dgvMerged     9.82e-15    -2.79e-37     Database of Genomic Variants: Structural Variant Regions (CNV, Inversion, In/del)
## 5  dgvSupporting     9.82e-15    -2.79e-37 Database of Genomic Variants: Supporting Structural Variants (CNV, Inversion, In/del)
## 6            gad     4.44e-20     4.94e-71                         Genetic Association Studies of Complex Diseases and Disorders
## 7        genscan     1.65e-12     1.10e-16                                                              Genscan Gene Predictions
## 8     kgProtMap2     9.87e-14     1.99e-59                                                                                  <NA>
## 9         rgdQtl     1.67e-23   -1.90e-296                                               Human Quantitative Trait Locus from RGD
## 10  ucscGenePfam     4.81e-18    -1.85e-68                                                            Pfam Domains in UCSC Genes
```

<img src="img/TAB123.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 26"
##               Row.names pval.TAB1.AA pval.TAB2.EA                                                                       V2
## 1           all_bacends    -5.51e-26    -2.51e-92                                                                     <NA>
## 2  altSeqLiftOverPslP10    -2.34e-05    -6.40e-25                           GRCh37 Alternate Sequence Lift Over Alignments
## 3   altSeqLiftOverPslP9    -2.34e-05    -6.40e-25                                                                     <NA>
## 4      altSeqPatchesP10    -3.58e-05    -1.75e-24                                     Patches to GRCh37 Reference Sequence
## 5            jaxQtlAsIs     7.53e-89    -8.00e-74                                  MGI Mouse QTLs Coarsely Mapped to Human
## 6           laminB1Lads    -4.14e-04    -1.45e-23                         NKI LADs (Lamina Associated Domains, Tig3 cells)
## 7          orfeomeGenes    -2.93e-18   -2.24e-100                                                                     <NA>
## 8           orfeomeMrna    -2.93e-18   -2.24e-100                                        ORFeome Collaboration Gene Clones
## 9              pubsBlat    -2.21e-35   -2.31e-106                        Sequences in Articles: PubmedCentral and Elsevier
## 10          pubsBlatPsl    -3.49e-23    -3.79e-59 Individual Sequence Matches of One Selected Article from Sequences Track
```

<img src="img/TAB124.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-altSplicing analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly DEPLETED associations is: 1"
##       Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1 strangeSplice     2.94e-06    -1.00e+00 <NA>
```

<img src="img/TAB125.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-chromStates analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##           Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1 10_Txn_Elongation     3.16e-01    -6.05e-14 <NA>
```

<img src="img/TAB126.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 4"
##          Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1      11_Weak_Txn     1.42e-03    -1.87e-01 <NA>
## 2 13_Heterochromlo    -4.37e-05     2.37e-52 <NA>
## 3  7_Weak_Enhancer    -1.79e-01     1.19e-09 <NA>
## 4 9_Txn_Transition    -6.37e-03    -8.76e-01 <NA>
```

<img src="img/TAB127.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-coriellCNVs analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##    Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1 Fibroblast    -2.07e-83     4.28e-47 <NA>
```

<img src="img/TAB128.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 2"
##                    Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1               B_Lymphocyte    -2.43e-25   -2.69e-143 <NA>
## 2 Chorionic_villus_cell_line    3.32e-130    -6.02e-31 <NA>
```

<img src="img/TAB129.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-gapLocations analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-genomicVariants analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly DEPLETED associations is: 6"
##     Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1     Complex    -2.66e-01    -8.12e-07 <NA>
## 2    Deletion    -1.39e-04    -6.22e-08 <NA>
## 3 Duplication    -2.72e-04    -3.96e-22 <NA>
## 4        Gain    -1.13e-25   -6.68e-110 <NA>
## 5   Inversion    -5.26e-12    -2.15e-56 <NA>
## 6        Loss     9.45e-31    -8.12e-07 <NA>
```

<img src="img/TAB1210.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-H3K4me3 analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 21"
##                          Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1           H3K4me3_Adipose_Nuclei     1.32e-02    -8.89e-10 <NA>
## 2             H3K4me3_Adult_Kidney     1.32e-02    -2.86e-05 <NA>
## 3       H3K4me3_CD19_Primary_Cells     1.42e-02    -1.32e-04 <NA>
## 4      H3K4me3_CD34_Cultured_Cells     5.18e-03    -5.28e-05 <NA>
## 5  H3K4me3_CD4_Naive_Primary_Cells     2.02e-03    -1.63e-04 <NA>
## 6  H3K4me3_CD8_Naive_Primary_Cells     2.83e-05    -2.54e-05 <NA>
## 7      H3K4me3_Colon_Smooth_Muscle     1.03e-03    -5.05e-08 <NA>
## 8           H3K4me3_Colonic_Mucosa     2.89e-02    -3.12e-05 <NA>
## 9     H3K4me3_Rectal_Smooth_Muscle     6.77e-05    -8.29e-09 <NA>
## 10          H3K4me3_Stomach_Mucosa     4.48e-03    -2.22e-04 <NA>
```

<img src="img/TAB1211.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 10"
##                               Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1        H3K4me3_Brain_Anterior_Caudate     7.79e-13    -3.14e-03 <NA>
## 2         H3K4me3_Brain_Cingulate_Gyrus     2.65e-07    -1.58e-03 <NA>
## 3      H3K4me3_Brain_Hippocampus_Middle     1.38e-09    -2.85e-03 <NA>
## 4  H3K4me3_Brain_Inferior_Temporal_Lobe     1.06e-12    -2.86e-05 <NA>
## 5        H3K4me3_Brain_Mid_Frontal_Lobe     9.33e-08    -1.32e-04 <NA>
## 6        H3K4me3_Brain_Substantia_Nigra     3.63e-08    -1.29e-03 <NA>
## 7             H3K4me3_CD3_Primary_Cells     1.64e-17    -1.98e-03 <NA>
## 8            H3K4me3_CD34_Primary_Cells     1.88e-06    -3.11e-03 <NA>
## 9               H3K4me3_Skeletal_Muscle     2.20e-20    -1.67e-07 <NA>
## 10           H3K4me3_Treg_Primary_Cells    -2.34e-02     3.00e-08 <NA>
```

<img src="img/TAB1212.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-ncRnas analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-repeats analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly DEPLETED associations is: 2"
##   Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1      LINE    -5.19e-02     4.67e-07 <NA>
## 2      SINE     5.76e-14    -1.00e+00 <NA>
```

<img src="img/TAB1213.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-tfbsConserved analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-tfbsEncode analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 3"
##   Row.names pval.TAB1.AA pval.TAB2.EA   V2
## 1      CTCF     5.59e-01    -1.30e-03 <NA>
## 2    POLR2A     1.00e+00    -7.88e-09 <NA>
## 3       YY1     9.96e-01    -6.02e-04 <NA>
```

<img src="img/TAB1214.png" title="plot of chunk TAB12" alt="plot of chunk TAB12" width="700" />

TAB1-AA analysis
================


```
## [1] "---------------------------------------------------------------"
## [1] "1Encode analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 369"
##                                              Row.names   TAB1.AA                                                                     V2
## 1              wgEncodeBroadHistoneH1hescH3k27me3StdPk 2.48e-165      H1-hESC H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2               wgEncodeBroadHistoneH1hescH3k4me1StdPk 1.32e-114       H1-hESC H3K4me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3               wgEncodeBroadHistoneH1hescSap3039731Pk 1.31e-135 H1-hESC SAP30 (39731) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4                wgEncodeBroadHistoneHsmmH3k27me3StdPk 1.85e-145         HSMM H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5                 wgEncodeBroadHistoneNhaH3k27me3StdPk 4.24e-143         NH-A H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6                  wgEncodeBroadHistoneNhdfadCtcfStdPk 1.64e-167          NHDF-Ad CTCF Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 7    wgEncodeSunyAlbanyGeneStGm12878Igf2bp1RbpAssocRna 4.95e-153                                                                   <NA>
## 8  wgEncodeSunyAlbanyGeneStGm12878Igf2bp1RbpAssocRnaV2 4.95e-153 GM12878 IGF2BP1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## 9         wgEncodeSunyAlbanyGeneStK562Celf1RbpAssocRna 2.92e-128                                                                   <NA>
## 10      wgEncodeSunyAlbanyGeneStK562Celf1RbpAssocRnaV2 2.92e-128      K562 CELF1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB11.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 132"
##                                              Row.names    TAB1.AA                                                                       V2
## 1               wgEncodeBroadHistoneK562Ezh239875StdPk  -9.47e-57       K562 EZH2 (39875) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2                 wgEncodeBroadHistoneK562H3k9me1StdPk  -1.09e-77            K562 H3K9me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3                wgEncodeBroadHistoneK562H4k20me1StdPk  -1.74e-57           K562 H4K20me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4                       wgEncodeBroadHistoneK562NcorPk  -4.50e-53               K562 NCoR Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5                    wgEncodeBroadHistoneK562P300StdPk  -8.70e-59               K562 P300 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6                       wgEncodeBroadHistoneK562PcafPk  -8.70e-59               K562 PCAF Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 7           wgEncodeGisChiaPetMcf7CtcfInteractionsRep2 -9.06e-172              MCF-7 CTCF ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 8           wgEncodeGisChiaPetMcf7Pol2InteractionsRep2  -2.42e-78              MCF-7 Pol2 ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 9    wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRna -3.87e-102                                                                     <NA>
## 10 wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRnaV2 -3.87e-102 HeLa-S3 RIP-Input RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB12.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "2other analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 42"
##        Row.names  TAB1.AA                                                                                    V2
## 1        affyU95 8.67e-16                              Alignments of Affymetrix Consensus/Exemplars from HG-U95
## 2       ccdsGene 1.12e-24                                                                         Consensus CDS
## 3      dgvMerged 9.82e-15     Database of Genomic Variants: Structural Variant Regions (CNV, Inversion, In/del)
## 4  dgvSupporting 9.82e-15 Database of Genomic Variants: Supporting Structural Variants (CNV, Inversion, In/del)
## 5     fishClones 1.06e-25                                           Clones Placed on Cytogenetic Map Using FISH
## 6            gad 4.44e-20                         Genetic Association Studies of Complex Diseases and Disorders
## 7     jaxQtlAsIs 7.53e-89                                               MGI Mouse QTLs Coarsely Mapped to Human
## 8         rgdQtl 1.67e-23                                               Human Quantitative Trait Locus from RGD
## 9   simpleRepeat 1.93e-14                                                          Simple Tandem Repeats by TRF
## 10  ucscGenePfam 4.81e-18                                                            Pfam Domains in UCSC Genes
```

<img src="img/TAB13.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 10"
##               Row.names   TAB1.AA                                                                       V2
## 1           all_bacends -5.51e-26                                                                     <NA>
## 2  altSeqLiftOverPslP10 -2.34e-05                           GRCh37 Alternate Sequence Lift Over Alignments
## 3   altSeqLiftOverPslP9 -2.34e-05                                                                     <NA>
## 4      altSeqPatchesP10 -3.58e-05                                     Patches to GRCh37 Reference Sequence
## 5           laminB1Lads -4.14e-04                         NKI LADs (Lamina Associated Domains, Tig3 cells)
## 6      lrgTranscriptAli -3.10e-03               Locus Reference Genomic (LRG) Fixed Transcript Annotations
## 7          orfeomeGenes -2.93e-18                                                                     <NA>
## 8           orfeomeMrna -2.93e-18                                        ORFeome Collaboration Gene Clones
## 9              pubsBlat -2.21e-35                        Sequences in Articles: PubmedCentral and Elsevier
## 10          pubsBlatPsl -3.49e-23 Individual Sequence Matches of One Selected Article from Sequences Track
```

<img src="img/TAB14.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "altSplicing analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##       Row.names  TAB1.AA   V2
## 1 strangeSplice 2.94e-06 <NA>
```

<img src="img/TAB15.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "chromStates analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##     Row.names  TAB1.AA   V2
## 1 11_Weak_Txn 1.42e-03 <NA>
```

<img src="img/TAB16.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 2"
##          Row.names   TAB1.AA   V2
## 1 13_Heterochromlo -4.37e-05 <NA>
## 2 9_Txn_Transition -6.37e-03 <NA>
```

<img src="img/TAB17.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "coriellCNVs analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##                    Row.names   TAB1.AA   V2
## 1 Chorionic_villus_cell_line 3.32e-130 <NA>
```

<img src="img/TAB18.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 2"
##      Row.names   TAB1.AA   V2
## 1 B_Lymphocyte -2.43e-25 <NA>
## 2   Fibroblast -2.07e-83 <NA>
```

<img src="img/TAB19.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "gapLocations analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "genomicVariants analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##   Row.names  TAB1.AA   V2
## 1      Loss 9.45e-31 <NA>
```

<img src="img/TAB110.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 4"
##     Row.names   TAB1.AA   V2
## 1    Deletion -1.39e-04 <NA>
## 2 Duplication -2.72e-04 <NA>
## 3        Gain -1.13e-25 <NA>
## 4   Inversion -5.26e-12 <NA>
```

<img src="img/TAB111.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "H3K4me3 analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 15"
##                               Row.names  TAB1.AA   V2
## 1        H3K4me3_Brain_Anterior_Caudate 7.79e-13 <NA>
## 2         H3K4me3_Brain_Cingulate_Gyrus 2.65e-07 <NA>
## 3      H3K4me3_Brain_Hippocampus_Middle 1.38e-09 <NA>
## 4  H3K4me3_Brain_Inferior_Temporal_Lobe 1.06e-12 <NA>
## 5        H3K4me3_Brain_Mid_Frontal_Lobe 9.33e-08 <NA>
## 6        H3K4me3_Brain_Substantia_Nigra 3.63e-08 <NA>
## 7             H3K4me3_CD3_Primary_Cells 1.64e-17 <NA>
## 8            H3K4me3_CD34_Primary_Cells 1.88e-06 <NA>
## 9       H3K4me3_CD8_Naive_Primary_Cells 2.83e-05 <NA>
## 10              H3K4me3_Skeletal_Muscle 2.20e-20 <NA>
```

<img src="img/TAB112.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "ncRnas analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "repeats analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##   Row.names  TAB1.AA   V2
## 1      SINE 5.76e-14 <NA>
```

<img src="img/TAB113.png" title="plot of chunk TAB1" alt="plot of chunk TAB1" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "tfbsConserved analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "tfbsEncode analysis"
## [1] "---------------------------------------------------------------"
```

TAB1-EA analysis
================


```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-1Encode analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 369"
##                                              Row.names   TAB1.AA                                                                     V2
## 1              wgEncodeBroadHistoneH1hescH3k27me3StdPk 2.48e-165      H1-hESC H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2               wgEncodeBroadHistoneH1hescH3k4me1StdPk 1.32e-114       H1-hESC H3K4me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3               wgEncodeBroadHistoneH1hescSap3039731Pk 1.31e-135 H1-hESC SAP30 (39731) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4                wgEncodeBroadHistoneHsmmH3k27me3StdPk 1.85e-145         HSMM H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5                 wgEncodeBroadHistoneNhaH3k27me3StdPk 4.24e-143         NH-A H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6                  wgEncodeBroadHistoneNhdfadCtcfStdPk 1.64e-167          NHDF-Ad CTCF Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 7    wgEncodeSunyAlbanyGeneStGm12878Igf2bp1RbpAssocRna 4.95e-153                                                                   <NA>
## 8  wgEncodeSunyAlbanyGeneStGm12878Igf2bp1RbpAssocRnaV2 4.95e-153 GM12878 IGF2BP1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## 9         wgEncodeSunyAlbanyGeneStK562Celf1RbpAssocRna 2.92e-128                                                                   <NA>
## 10      wgEncodeSunyAlbanyGeneStK562Celf1RbpAssocRnaV2 2.92e-128      K562 CELF1 RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB21.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 132"
##                                              Row.names    TAB1.AA                                                                       V2
## 1               wgEncodeBroadHistoneK562Ezh239875StdPk  -9.47e-57       K562 EZH2 (39875) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2                 wgEncodeBroadHistoneK562H3k9me1StdPk  -1.09e-77            K562 H3K9me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3                wgEncodeBroadHistoneK562H4k20me1StdPk  -1.74e-57           K562 H4K20me1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4                       wgEncodeBroadHistoneK562NcorPk  -4.50e-53               K562 NCoR Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5                    wgEncodeBroadHistoneK562P300StdPk  -8.70e-59               K562 P300 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6                       wgEncodeBroadHistoneK562PcafPk  -8.70e-59               K562 PCAF Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 7           wgEncodeGisChiaPetMcf7CtcfInteractionsRep2 -9.06e-172              MCF-7 CTCF ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 8           wgEncodeGisChiaPetMcf7Pol2InteractionsRep2  -2.42e-78              MCF-7 Pol2 ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 9    wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRna -3.87e-102                                                                     <NA>
## 10 wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRnaV2 -3.87e-102 HeLa-S3 RIP-Input RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
```

<img src="img/TAB22.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-2other analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 42"
##        Row.names  TAB1.AA                                                                                    V2
## 1        affyU95 8.67e-16                              Alignments of Affymetrix Consensus/Exemplars from HG-U95
## 2       ccdsGene 1.12e-24                                                                         Consensus CDS
## 3      dgvMerged 9.82e-15     Database of Genomic Variants: Structural Variant Regions (CNV, Inversion, In/del)
## 4  dgvSupporting 9.82e-15 Database of Genomic Variants: Supporting Structural Variants (CNV, Inversion, In/del)
## 5     fishClones 1.06e-25                                           Clones Placed on Cytogenetic Map Using FISH
## 6            gad 4.44e-20                         Genetic Association Studies of Complex Diseases and Disorders
## 7     jaxQtlAsIs 7.53e-89                                               MGI Mouse QTLs Coarsely Mapped to Human
## 8         rgdQtl 1.67e-23                                               Human Quantitative Trait Locus from RGD
## 9   simpleRepeat 1.93e-14                                                          Simple Tandem Repeats by TRF
## 10  ucscGenePfam 4.81e-18                                                            Pfam Domains in UCSC Genes
```

<img src="img/TAB23.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 10"
##               Row.names   TAB1.AA                                                                       V2
## 1           all_bacends -5.51e-26                                                                     <NA>
## 2  altSeqLiftOverPslP10 -2.34e-05                           GRCh37 Alternate Sequence Lift Over Alignments
## 3   altSeqLiftOverPslP9 -2.34e-05                                                                     <NA>
## 4      altSeqPatchesP10 -3.58e-05                                     Patches to GRCh37 Reference Sequence
## 5           laminB1Lads -4.14e-04                         NKI LADs (Lamina Associated Domains, Tig3 cells)
## 6      lrgTranscriptAli -3.10e-03               Locus Reference Genomic (LRG) Fixed Transcript Annotations
## 7          orfeomeGenes -2.93e-18                                                                     <NA>
## 8           orfeomeMrna -2.93e-18                                        ORFeome Collaboration Gene Clones
## 9              pubsBlat -2.21e-35                        Sequences in Articles: PubmedCentral and Elsevier
## 10          pubsBlatPsl -3.49e-23 Individual Sequence Matches of One Selected Article from Sequences Track
```

<img src="img/TAB24.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-altSplicing analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##       Row.names  TAB1.AA   V2
## 1 strangeSplice 2.94e-06 <NA>
```

<img src="img/TAB25.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-chromStates analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##     Row.names  TAB1.AA   V2
## 1 11_Weak_Txn 1.42e-03 <NA>
```

<img src="img/TAB26.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 2"
##          Row.names   TAB1.AA   V2
## 1 13_Heterochromlo -4.37e-05 <NA>
## 2 9_Txn_Transition -6.37e-03 <NA>
```

<img src="img/TAB27.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-coriellCNVs analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##                    Row.names   TAB1.AA   V2
## 1 Chorionic_villus_cell_line 3.32e-130 <NA>
```

<img src="img/TAB28.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 2"
##      Row.names   TAB1.AA   V2
## 1 B_Lymphocyte -2.43e-25 <NA>
## 2   Fibroblast -2.07e-83 <NA>
```

<img src="img/TAB29.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-gapLocations analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-genomicVariants analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##   Row.names  TAB1.AA   V2
## 1      Loss 9.45e-31 <NA>
```

<img src="img/TAB210.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "The total number of significantly DEPLETED associations is: 4"
##     Row.names   TAB1.AA   V2
## 1    Deletion -1.39e-04 <NA>
## 2 Duplication -2.72e-04 <NA>
## 3        Gain -1.13e-25 <NA>
## 4   Inversion -5.26e-12 <NA>
```

<img src="img/TAB211.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-H3K4me3 analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 15"
##                               Row.names  TAB1.AA   V2
## 1        H3K4me3_Brain_Anterior_Caudate 7.79e-13 <NA>
## 2         H3K4me3_Brain_Cingulate_Gyrus 2.65e-07 <NA>
## 3      H3K4me3_Brain_Hippocampus_Middle 1.38e-09 <NA>
## 4  H3K4me3_Brain_Inferior_Temporal_Lobe 1.06e-12 <NA>
## 5        H3K4me3_Brain_Mid_Frontal_Lobe 9.33e-08 <NA>
## 6        H3K4me3_Brain_Substantia_Nigra 3.63e-08 <NA>
## 7             H3K4me3_CD3_Primary_Cells 1.64e-17 <NA>
## 8            H3K4me3_CD34_Primary_Cells 1.88e-06 <NA>
## 9       H3K4me3_CD8_Naive_Primary_Cells 2.83e-05 <NA>
## 10              H3K4me3_Skeletal_Muscle 2.20e-20 <NA>
```

<img src="img/TAB212.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-ncRnas analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-repeats analysis"
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly ENRICHED associations is: 1"
##   Row.names  TAB1.AA   V2
## 1      SINE 5.76e-14 <NA>
```

<img src="img/TAB213.png" title="plot of chunk TAB2" alt="plot of chunk TAB2" width="700" />

```
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-tfbsConserved analysis"
## [1] "---------------------------------------------------------------"
## [1] "---------------------------------------------------------------"
## [1] "data/TAB1-AA-tfbsEncode analysis"
## [1] "---------------------------------------------------------------"
```
