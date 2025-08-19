# Leptospira LIC13259 C8G screen
 
This repo documents the code used to carry out the analysis for the paper titled: 'AlphaFold reveals how pathogenic Leptospira use cross-kingdom thiol-disulphide exchange to evade the complement membrane attack complex'.

<center>
![image](https://github.com/CBFLivUni/Leptospira-LIC13259-C8G/blob/main/summary_figure.png)
</center>

## Set-up 

Analysis was carried out on the Barkla1 cluster. It should be fully reproducible starting from two fasta files: 'LIC13259_sequences_no_his_tag.fasta' and 'c8_uniprotkb_2024_07_25.fasta' contained within the `raw_fasta` folder.

In addition to this, the analysis requires three key pieces of software to be installed (in addition to a standard Rstudio set-up): 
- [ColabFold](https://github.com/sokrypton/ColabFold) installation instructions can be found here
- [FoldSeek](https://github.com/steineggerlab/foldseek) installation instructions and tutorials can be found here
- [Alphapulldown](https://github.com/KosinskiLab/AlphaPulldown) installation instructions and tutorials

## Running the analysis 

The analysis was carried out using a series of bash scripts and quarto markdown documents (qmds). These should be run in order. Supporting scripts are found inside the `scripts` folder.
