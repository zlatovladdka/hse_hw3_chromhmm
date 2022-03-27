# hse_hw3_chromhmm
## Colab
https://colab.research.google.com/drive/1dVrTvP_78V-dzHprN3I5Ef9H2ZrgZE18?usp=sharing

```bash
!curl -O https://raw.githubusercontent.com/deepjavalibrary/d2l-java/master/tools/fix-colab-gpu.sh && bash fix-colab-gpu.sh
!curl -O https://raw.githubusercontent.com/deepjavalibrary/d2l-java/master/tools/colab_build.sh && bash colab_build.sh
!java --list-modules | grep "jdk.jshell"
!wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
!unzip /content/ChromHMM.zip
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562ControlStdAlnRep1.bam -O Control.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562Suz12051317AlnRep1.bam   -O Suz.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k36me3StdAlnRep1.bam  -O H3k36.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k9me1StdAlnRep1.bam -O H3k9.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562CtcfStdAlnRep2.bam -O Ctcf.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562Cbx2AlnRep1.bam -O Cbx.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H2azStdAlnRep2.bam -O H2az.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562Pol2bStdAlnRep1.bam -O Pol.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k79me2StdAlnRep2.bam  -O H3k79.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k4me1StdAlnRep1.bam -O H3k4.bam
!wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k27me3StdAlnRep1.bam -O H3k27.bam
!java -mx5000M -jar /content/ChromHMM/ChromHMM.jar BinarizeBam -b 200  /content/ChromHMM/CHROMSIZES/hg19.txt /content/bam_files cellmarkfiletable.txt   binarizedData
!java -mx1600m -jar /content/ChromHMM/ChromHMM.jar LearnModel -p 0 binarizedData output_data 10 hg19
from google.colab import drive
drive.mount('/content/gdrive')
!mv /content/output_data /content/gdrive/MyDrive/HW3_chromhmm
```

## Данные
Клеточная линия K562

Гистоновая метка | Имя файла | Используемый файл
--- | --- | ---
Control | wgEncodeBroadHistoneK562ControlStdAlnRep1.bam | Control.bam
Suz12 | wgEncodeBroadHistoneK562Suz12051317AlnRep1.bam | Suz.bam
H3k36me3 | wgEncodeBroadHistoneK562H3k36me3StdAlnRep1.bam | H3k36.bam
H3k9me1 | wgEncodeBroadHistoneK562H3k9me1StdAlnRep1.bam | H3k9.bam
CTCF | wgEncodeBroadHistoneK562CtcfStdAlnRep2.bam | Ctcf.bam
CBX2 | wgEncodeBroadHistoneK562Cbx2AlnRep1.bam | Cbx.bam
H2A.Z | wgEncodeBroadHistoneK562H2azStdAlnRep2.bam | H2az.bam
Pol2B | wgEncodeBroadHistoneK562Pol2bStdAlnRep1.bam | Pol.bam
H3K79me2 | wgEncodeBroadHistoneK562H3k79me2StdAlnRep2.bam | H3K79.bam
H3K4me1 | wgEncodeBroadHistoneK562H3k4me1StdAlnRep1.bam | H3K4.bam
H3K27me3 | wgEncodeBroadHistoneK562H3k27me3StdAlnRep1.bam | H3K27.bam


## Отчет CromHMM

Fold Enrichment overlap | Fold Enrichment TES neighborhood | Fold Enrichment TSS neighborhood | Emission Parameters | Transition parameters
--- | --- | --- | --- | ---
![](/output/K562_10_overlap.png) |  ![](/output/K562_10_RefSeqTES_neighborhood.png) | ![](/output/K562_10_RefSeqTSS_neighborhood.png) |  ![](/output/emissions_10.png) | ![](/output/transitions_10.png)


## UCSC 
Визуализация данных в геномном браузере.

![](/output/ucsc1.png)
![](/output/ucsc2.png)


## Бонус

```python
import pandas as pd
df = pd.read_csv('/content/gdrive/MyDrive/HW3_chromhmm/output_data/K562_10_dense.bed', encoding='utf-8', sep='\t', comment='t', header=None)
df.loc[df[3] == 1, 3] = '1_insulator'
df.loc[df[3] == 2, 3] = '2_promoter'
df.loc[df[3 ]== 3, 3] = '3_poised_promoter'
df.loc[df[3] == 4, 3] = '4_rpr_heterochromatin'
df.loc[df[3] == 5, 3] = '5_trx_elongation'
df.loc[df[3] == 6, 3] = '6_skip_exons'
df.loc[df[3] == 7, 3] = '7_skip_exons'
df.loc[df[3] == 8, 3] = '8_trx_elongation'
df.loc[df[3] == 9, 3] = '9_promoter'
df.loc[df[3] == 10, 3] = '10_promoter'
df.to_csv('/content/gdrive/MyDrive/HW3_chromhmm/K562_10_dense_edit.bed', sep='\t', index=False, header=None)
```


