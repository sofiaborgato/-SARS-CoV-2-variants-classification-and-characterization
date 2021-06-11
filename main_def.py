#!/usr/bin/env python
# coding: utf-8

import os
path=os.getcwd() 

print('Downloading the Git repo...')
import git
git.Git(path).clone('https://github.com/sofiaborgato/Bioinformatic_project.git')

os.chdir(os.path.join(path,'Bioinformatic_project'))

print(os.getcwd())

import sys

sys.path.append(os.path.join(path,'Bioinformatic_project'))
from preprocessing_utils import *
from supervised_classification_utils  import *
from clustering_utils import nullscan
from clustering_utils import clustering
import tarfile

# open file
file = tarfile.open('./data.tar.xz') 
# extracting file
file.extractall()
  
file.close()

# open file
file = tarfile.open('./data_local_concat.tar.xz') 
# extracting file
file.extractall()
  
file.close()


# In[4]:


import pandas as pd
total_al = pd.read_csv('./data_local_concat/local_aligned_.csv')
total_stats=   pd.read_csv('./data_local_concat/stats.csv')
genomes_aligned_class=pd.read_csv('./demo_class_aligned_new.csv')
stats_class=pd.read_csv('./demo_class_stats.csv')
key_mut_clust=pd.read_csv('./demo_class_key_mutations.csv')

genomes_aligned_clust=pd.read_csv('./demo_clust_aligned_new.csv')
stats_clust=pd.read_csv('./demo_clust_stats.csv')
key_mut_clust=pd.read_csv('./demo_clust_key_mutations.csv')

# In[5]:


in_path = input("Please enter the path of your fasta file, otherwise type 0 to use demo data --> ")
while in_path[-6:]!='.fasta' and in_path!='0':
    print("The format of the file isn't correct. Please insert a .fasta file")
    in_path = input("Please enter the path of your fasta file, otherwise type 0 to use demo data --> ")


an_type = input("\nWhich kind of analysis do you desire to perform? \n 1. Supervised classification : Labels your samples according to the known variants. More accurate and less sensitive to the outliers (Advised for files with <100 samples); \n 2. Unsupervised clustering : More sensitive to the outliers but capable of highlighting new lineages of the genome (Advised for files with >100 samples). \n\n Insert your choice (1/2) ---> ")
while an_type != '1' and an_type != '2':
    print("TPlease type 1 or 2.")
    an_type = input("\nWhich kind of analysis do you desire to perform? \n 1. Supervised classification : Labels your samples according to the known variants. More accurate and less sensitive to the outliers (Advised for files with <100 samples); \n 2. Unsupervised clustering : More sensitive to the outliers but capable of highlighting new lineages of the genome (Advised for files with >100 samples). \n\n Insert your choice (1/2) ---> ")



if in_path != '0':
    stats, genomes_aligned, mutation_list = align_and_process(in_path,string_length=500)
    key_mut = key_mutations(mutation_list, len(stats))
    export_data(stats, genomes_aligned,key_mut, name = "new")

if in_path ==  '0' and an_type == '1':
    stats = stats_class
    genomes_aligned = genomes_aligned_class
    key_mut = key_mut_class
    export_data(stats, genomes_aligned,key_mut, name = "./Output/demo")
    prediction=classifier(total_stats,stats)
    
    
    
elif in_path == '0' and an_type == '2':
	stats = stats_clust
	genomes_aligned = genomes_aligned_clust
	key_mut = key_mut_clust
	export_data(stats, genomes_aligned,key_mut, name = "./Output/demo")
	prediction, num_new_var = clustering(total_stats,stats)
	stats['Predicted'] = prediction
	if num_new_var > 0:
		for i in range(num_new_var):
			test_genomes = genomes_aligned[stats.Predicted == i+6]
			mutation_list_new, _ = process_mutations(test_genomes['Reference aligned'], test_genomes['Full sequence'])
			key_mut_new = key_mutations(mutation_list_new, len(test_genomes))
			name = "New variant " + str(i+1) + "key mutations.csv"
			key_mut_new.to_csv(name)

