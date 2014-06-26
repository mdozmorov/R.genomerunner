








SLE-IC SNP sets analysis, ENCODE genome annotation data
========================================================
 

Double lines (===) mark the names of the regions being analyzed. 

Pluses (+++) mark the names of categories of genomic features, used for the enrichment analysis

Cell type-specific sets of epigenomic elements from the ENCODE project were used.

For each analysis, only the significant enrichment results (p < 0.01) are shown. If an analysis does not have any significant results, it is not displayed at all. The analyses are separated by single lines (---).

Negative p-value indicates depletion.





```
## [1] "==============================================================="
## [1] "SNP set analyzed: data.enc/TAB12-combined-ENC"
## [1] "==============================================================="
## [1] "The total number of significantly ENRICHED associations is: 223 Top 10, or less, are shown"
##                                     Row.names Enrich.pval                                                                     V2
## 1        wgEncodeBroadHistoneA549H2azEtoh02Pk   2.78e-148 A549 EtOH 0.02% H2A.Z Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2             wgEncodeBroadHistoneDnd41H2azPk   3.87e-174           Dnd41 H2A.Z Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 3      wgEncodeBroadHistoneHepg2H3k27me3StdPk   2.00e-207        HepG2 H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 4       wgEncodeBroadHistoneHsmmH3k27me3StdPk   3.47e-183         HSMM H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 5      wgEncodeBroadHistoneHuvecH3k27me3StdPk   2.25e-303        HUVEC H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 6              wgEncodeBroadHistoneK562RestPk   9.28e-156             K562 REST Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 7           wgEncodeBroadHistoneNhaH3k09me3Pk   2.86e-133          NH-A H3K9me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 8         wgEncodeBroadHistoneOsteoH3k27me3Pk   1.32e-298  Osteoblasts H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 9  wgEncodeGisChiaPetMcf7EraaInteractionsRep1   1.70e-211       MCF-7 ERalpha a ChIA-PET Interactions Rep 1 from ENCODE/GIS-Ruan
## 10     wgEncodeHaibGenotypeMyometrRegionsRep1   6.54e-157              Myometr Copy number variants Replicate 1 from ENCODE/HAIB
## [1] "---------------------------------------------------------------"
## [1] "The total number of significantly DEPLETED associations is: 389 Top 10, or less, are shown"
##                                             Row.names Enrich.pval                                                                       V2
## 1             wgEncodeBroadHistoneNhdfadH3k36me3StdPk  -1.91e-163        NHDF-Ad H3K36me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## 2       wgEncodeCshlLongRnaSeqA549CytosolPapJunctions  -6.50e-176            A549 cytosol polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## 3         wgEncodeCshlLongRnaSeqHuvecCellPapJunctions  -1.88e-202        HUVEC whole cell polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## 4       wgEncodeCshlLongRnaSeqK562NucleusPapJunctions  -2.76e-214            K562 nucleus polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## 5          wgEncodeGisChiaPetMcf7Pol2InteractionsRep2  -2.25e-303              MCF-7 Pol2 ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## 6           wgEncodeGisChiaPetNb4Pol2InteractionsRep1  -1.50e-206                NB4 Pol2 ChIA-PET Interactions Rep 1 from ENCODE/GIS-Ruan
## 7         wgEncodeGisRnaPetNhekNucleusPapClustersRep1  -1.95e-163    NHEK nucleus polyA+ clone-free RNA PET Clusters Rep 1 from ENCODE/GIS
## 8 wgEncodeSunyAlbanyGeneStHelas3RipinputRbpAssocRnaV2  -2.61e-267 HeLa-S3 RIP-Input RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## [1] "---------------------------------------------------------------"
```


```
## [1] "Top 20 up in condition 1"
##                                                     TAB1.AA TAB2.EA                                                                           V2
## wgEncodeBroadHistoneH1hescH3k27me3StdPk             164.606  -15.01            H1-hESC H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneHelas3Pol2bStdPk                 32.606 -225.71                HeLa-S3 Pol2 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeCshlLongRnaSeqA549CytosolPapJunctions         4.204 -212.22                A549 cytosol polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqBjCellPamJunctions              7.099 -208.00               BJ whole cell polyA- RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqHepg2CytosolPamJunctions        8.559 -171.33               HepG2 cytosol polyA- RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqHuvecCellPapJunctions           5.136 -240.89            HUVEC whole cell polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqK562NucleusPapJunctions         4.959 -252.24                K562 nucleus polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqNhekNucleusPamJunctions         5.147 -182.25                NHEK nucleus polyA- RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeCshlLongRnaSeqNhlfCellPapJunctions            5.698 -199.32             NHLF whole cell polyA+ RNA-seq Junctions Pooled from ENCODE/CSHL
## wgEncodeGisChiaPetMcf7Pol2InteractionsRep2          -77.615 -302.65                  MCF-7 Pol2 ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## wgEncodeGisRnaPetA549CellPapClusters                  7.304 -192.85    A549 whole cell polyA+ clone-free RNA PET Clusters Pooled from ENCODE/GIS
## wgEncodeGisRnaPetA549CytosolPapClusters               7.507 -194.35       A549 cytosol polyA+ clone-free RNA PET Clusters Pooled from ENCODE/GIS
## wgEncodeGisRnaPetA549NucleusPapClusters               8.352 -174.75       A549 nucleus polyA+ clone-free RNA PET Clusters Pooled from ENCODE/GIS
## wgEncodeGisRnaPetHepg2CytosolPapClustersRep1          5.129 -208.01      HepG2 cytosol polyA+ clone-based RNA PET Clusters Rep 1 from ENCODE/GIS
## wgEncodeGisRnaPetHepg2NucleusPapClustersRep1          0.000 -191.10      HepG2 nucleus polyA+ clone-based RNA PET Clusters Rep 1 from ENCODE/GIS
## wgEncodeGisRnaPetNhekNucleusPapClustersRep1           0.000 -214.25        NHEK nucleus polyA+ clone-free RNA PET Clusters Rep 1 from ENCODE/GIS
## wgEncodeGisRnaPetSknshCellPapClusters                 7.107 -172.29 SK-N-SH whole cell polyA+ clone-free RNA PET Clusters Pooled from ENCODE/GIS
## wgEncodeSunyAlbanyGeneStH1hescRipinputRbpAssocRnaV2   7.100 -207.30     H1-hESC RIP-Input RBP Associated RNA by RIP-chip GeneST from ENCODE/SUNY
## [1] "---------------------------------------------------------------"
## [1] "Top 20 up in condition 2"
##                                              TAB1.AA TAB2.EA                                                                      V2
## wgEncodeBroadHistoneA549H2azEtoh02Pk         -2.3934  206.92  A549 EtOH 0.02% H2A.Z Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneDnd41H2azPk              -3.3346  245.11            Dnd41 H2A.Z Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneHepg2H3k27me3StdPk        0.0000  271.33         HepG2 H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneHsmmH3k9me3StdPk        -11.0139  177.45           HSMM H3K9me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562Cbpsc369Pk          -42.4797  158.40      K562 CBP (sc-369) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562Cbx3sc101004Pk      -33.8696  194.22  K562 CBX3 (SC-101004) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562Ezh239875StdPk      -56.0238  112.20      K562 EZH2 (39875) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562H3k27me3StdPk        -7.7078  199.42          K562 H3K27me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562H3k9me3StdPk        -28.2037  207.55           K562 H3K9me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562NcorPk              -52.3472  145.51              K562 NCoR Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562P300StdPk           -58.0606  112.12              K562 P300 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562RestPk              -22.1835  271.91              K562 REST Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneK562Setdb1Pk            -41.4554  138.23            K562 SETDB1 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneMonocd14ro1746H2azPk     -0.7618  166.40  Monocytes CD14+ H2A.Z Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneNhaEzh239875Pk          -19.4592  158.37      NH-A EZH2 (39875) Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneNhaH3k09me3Pk           -10.3243  192.79           NH-A H3K9me3 Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeBroadHistoneOsteoP300kat3bPk        -30.3856  205.94 Osteoblasts P300 KAT3B Histone Mods by ChIP-seq Peaks from ENCODE/Broad
## wgEncodeGisChiaPetMcf7CtcfInteractionsRep2 -171.0430   13.28             MCF-7 CTCF ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## wgEncodeGisChiaPetMcf7EraaInteractionsRep1  -17.5848  302.65        MCF-7 ERalpha a ChIA-PET Interactions Rep 1 from ENCODE/GIS-Ruan
## wgEncodeGisChiaPetMcf7EraaInteractionsRep2  -35.0549  190.40        MCF-7 ERalpha a ChIA-PET Interactions Rep 2 from ENCODE/GIS-Ruan
## [1] "---------------------------------------------------------------"
```

<img src="img/TAB1vsTAB21.png" title="plot of chunk TAB1vsTAB2" alt="plot of chunk TAB1vsTAB2" width="700" /><img src="img/TAB1vsTAB22.png" title="plot of chunk TAB1vsTAB2" alt="plot of chunk TAB1vsTAB2" width="700" />

