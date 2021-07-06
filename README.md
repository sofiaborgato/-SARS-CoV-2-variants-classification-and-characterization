# SARS-CoV-2 variants classification and characterization
Preprocessing and description of SARS-CoV-2 genomes and classification into variants according to their key mutations with respect to the reference genome sequence Wuhan-Hu-1.
Final project for the course "Bioinformatics", A.Y. 2020/2021.

## Data
We used genome sequences from [`GISAID Dataset`](https://www.gisaid.org/). The dataset contains sequences 29kb long from the genome of the virus with <1% of unknown bases, divided according to the lineage considered and the date and place of submission. We downloaded 6500 sequences:
* 1000 samples from the original cases from Wuhan in January 2020;
* 1000 samples from lineage B.1.1.7 (English variant);
* 1000 samples from lineage B.1.351 (South-African variant);
* 1000 samples from lineage P.1 (Brazilian variant);
* 1000 samples from lineage B.1.427/B.1.429 (Californian variant);
* 1000 samples from lineage B.1.525 (Nigerian variant);
* 500 samples from lineage B.1.617 (Indian variant).
### Demo Data
For the demo trial we prepared two different dataset one for the supevised task and one for the clustering task.
The demo dataset for the classification task is composed by:

* 30 samples from the original cases from Wuhan in January 2020;
* 30 samples from lineage B.1.1.7 (English variant);
* 30 samples from lineage B.1.351 (South-African variant);
* 30 samples from lineage P.1 (Brazilian variant);
* 30 samples from lineage B.1.427/B.1.429 (Californian variant);
* 30 samples from lineage B.1.525 (Nigerian variant);
 
The demo dataset for the unsupervised task is composed by :

* 10 samples from lineage B.1.1.7 (English variant);
* 10 samples from lineage B.1.351 (South-African variant);
* 10 samples from lineage P.1 (Brazilian variant);
* 10 samples from lineage B.1.427/B.1.429 (Californian variant);
* 10 samples from lineage B.1.525 (Nigerian variant);
* 150 samples from lineage B.1.617 (Indian variant).



## Report of the work
A detailed description of the project and the analyses performed can be found [`here`](./Report/Report.pdf)

## Description of the repository
The root folder of the repo contains:
* the code to be run as explained in the **Instructions** section;
* the folder [`Data`](./Data), containing the FASTA files with the genome sequences and the inputs for the demo trial;
* the folder [`Report`](./Report) with a detailed description of the project and the analyses performed;
* the folder [`Output`](./Output) where the output of the pipeline is stored, which is then moved to the preferred location
* the folder [`Classifiers`](./Classifiers) containing the code for the different classifiers implemented;

## Instructions 
In order to run the code it's sufficient to download and run from terminal the file [`main_def.py`](main_def.py).
By following the instructions shown in the terminal you should be able to run two kinds of analysis:
 1. Supervised classification : Labels your samples according to the known variants. More accurate and less sensitive to the outliers (Advised for files with <100 samples);
 2. Unsupervised clustering : More sensitive to the outliers but capable of highlighting new lineages of the genome (Advised for files with >100 samples).
 
If you wish to analyse a FASTA file with new genomes sequences, it will be asked to specify the path to the file. Otherwise it's still possible to run the code with a demo dataset already aligned.
The outputs of the execution are stored in a specified output path.
 
 
