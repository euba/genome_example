﻿# Purpose
This setup and R script is to create a customizeable visualization of genome information.

## Genome annotation
As a prerequisite, you have to annotate your genome sequence with  [kbase](https://kbase.us/) with Annotation Apps such as RAST or Prokka. When done, export the annotation as json file on the export menu of your genome annotation.

## Categorize KO ids to your genes
Next, run the first part of the Rscript *circos_genome.R* to create a protein fasta file with the amino acid sequences for each gene. You can then upload this fasta file on [BlastKOALA](https://www.kegg.jp/blastkoala/) to annotate the KEGG Orthology (KO) ID for each gene. The output of this process (will take a while to send to you by the web service) is a text file with the gene id and the KO id.

## Formating and plotting
Now you basically have all of the things neccessary to plot your genome. Let me just quickly summarize again what you need:
* Genome sequence as fasta file
* Genome annotation as json file from [kbase](https://kbase.us/) 
* Gene to KO id translation table from [BlastKOALA](https://www.kegg.jp/blastkoala/) 

If you have all these items, you can proceed with the other parts of the Rscript *circos_genome.R*.
