# -*- coding: utf-8 -*-
"""preprocessing_utils.ipynb

Automatically generated by me.

Original file is located at
    https://colab.research.google.com/drive/1-TnDPz1h9PhMPnbH6HgXvz1po8ytP4Ej
"""



# Import pairwise2 module
#pip3 install Bio
from Bio import pairwise2
import os
import numpy as np
import pandas as pd
from Bio.Seq import Seq




#Function to replace bases different from A,C,G,T with X
def wrong_basis(text):
    chars = ['N','R','Y','K','M','S','W','B','D','H','V']
    for c in chars:
        text = text.replace(c,'X')
    return text

#Function to read and import a fasta file 
def read(path):
    reads=[] #array of the reads to be returned
    bases = ['A','C','G','T']
    first = True
    with open(path) as f:
        genome = ""
        for line in f:
            if line[0] != '>': #if the sample is split through multiple rows, the rows are concatenated into one single string until the next header row
                row=line.replace('\n','')
                genome = genome + row
            elif first:
                first = False
            else: #if the header row isn't the first of the fasta, the sample is stored in the list reads after being preprocessed
                genome=genome.upper()
                genome=wrong_basis(genome) 
                reads.append(genome)
                genome = ""

        genome=genome.upper()
        genome=wrong_basis(genome)         
        reads.append(genome)            
    
    return reads #the list of samples is returned

#Function to align two sequences with local alignment. 
#Due to high computational cost of working with very long rows, the task is divided in steps which are then concatenated properly
def align(X,Y, corr = 2, mis = -0.1, gap = -2, rgap = -0.2, span = 2000,step = 1000):
  span = span                 #the length of the partial strings to be aligned
  step = step                 #the length of the step at every iteration. It works best when is ~ half of the step
  up = int(span/2 + step/2)   #the upper index of the strings to be concatenated after being aligned
  low = int(span/2 - step/2)  #the lower index of the strings to be concatenated after being aligned

  #Scores
  corr_score = corr     #score for a correct match
  mis_score = mis   #score for a mismatch
  gap_score = gap   #score for a gap
  rgap_score = rgap #score for multiple gaps

  first = True
  seq1 = ""
  seq2 = ""

  counter = 0
  while ((counter*step + span) < len(Y)):

    #local alignment of the partial strings with a function from BioPython
    alignments = pairwise2.align.localms( X[step*counter: step*counter + span],  Y[step*counter: step*counter + span],corr_score,mis_score, gap_score, rgap_score, one_alignment_only=True)
    
    if first:
      #if this is the first partial alignment we consider the substrings from start to up index
      seq1 = seq1 + alignments[0][0][:up]
      seq2 = seq2 + alignments[0][1][:up]
      first = False
    else:
      #if this is a middle partial alignment we consider the substrings from low to up index
      seq1 = seq1 + alignments[0][0][low:up]
      seq2 = seq2 + alignments[0][1][low:up]

    counter = counter + 1
  #if this is the last partial alignment we consider the substrings from low index to the end
  alignments = pairwise2.align.localms( X[step*counter: ],  Y[step*counter: ],corr_score,mis_score, gap_score, rgap_score,one_alignment_only=True)
  seq1 = seq1 + alignments[0][0][low:]
  seq2 = seq2 + alignments[0][1][low:]

  #seq_s can be used to visualize the alignment correctly formatted if needed
  seq_s = ""
  c = len(seq1)
  for i in range(c):
    if seq1[i] == seq2[i]:
      seq_s = seq_s + "|"
    if seq1[i] == "-" or seq2[i] == "-":
      seq_s = seq_s + " "
    else:
      if seq1[i] != seq2[i]:
        seq_s = seq_s + "."

  return seq1,seq2,seq_s

#Function to create a dictionary from a list with its items and their frequencies
def listToFreq(list_m):
  freq = [list_m.count(p) for p in list_m]
  return dict(list(zip(list_m,freq)))

#Function to sort a dictionary according to the value
def sortFreqDict(freqdict):
  aux = [(freqdict[key], key) for key in freqdict]
  aux.sort()
  aux.reverse()
  return aux





def process_mutations(ref,seq_2):

  mutation_types = ['silent','non sense','missense conservative','missense non conservative','deletion in frame','insertion in frame','deletion frame-shift', 'insertion frame-shift']
  
  #Dictionaries to classify the amminoacids
  gen_code = {'TTT' : 'F' , 'TTC' : 'F' , 'TTA' : 'L' , 'TTG' : 'L' , 'CTT' : 'L' , 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L', 'ATT' : 'I', 'ATC' : 'I' , 'ATA' : 'I', 'ATG' : 'M', 'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V', 'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S', 'TCG' : 'S', 'CCT' : 'P', 'CCC' : 'P' , 'CCA' : 'P', 'CCG' : 'P', 'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T', 'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A', 'TAT' : 'Y', 'TAC' : 'Y', 'TAA' : 'STOP', 'TAG' : 'STOP', 'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q', 'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG':'K', 'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E', 'TGT' : 'C', 'TGC' : 'C', 'TGA' : 'STOP', 'TGG' : 'W', 'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R' , 'CGG' : 'R', 'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R' , 'AGG' : 'R', 'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G'}
  amm_class = {'G':'aliphatic','A':'aliphatic','V':'aliphatic','L':'aliphatic','I':'aliphatic', 'S':'hydroxil','C':'hydroxil','T':'hydroxil','M':'hydroxil','P':'cyclic', 'F':'aromatic','Y':'aromatic','W':'aromatic', 'H':'basic','K':'basic','R':'basic','D':'acidic','E':'acidic','N':'acidic','Q':'acidic'}
  
  #Counters of the mutations in each gene
  ORF1ab = [0,0,0,0,0,0,0]
  S = [0,0,0,0,0,0,0]
  ORF3a = [0,0,0,0,0,0,0]
  E  = [0,0,0,0,0,0,0]
  M = [0,0,0,0,0,0,0]
  ORF6  = [0,0,0,0,0,0,0]
  ORF7a = [0,0,0,0,0,0,0]
  ORF7b = [0,0,0,0,0,0,0]
  ORF8 = [0,0,0,0,0,0,0]
  N  = [0,0,0,0,0,0,0]
  ORF10 = [0,0,0,0,0,0,0]
  NON_COD = [0,0,0,0,0,0,0]

  #List of the mutation names
  mutation_list = []

  #Starting base
  i = 265

  while i < np.min([29674,len(seq_2)]):
      
      #If a basis in ref is different from the corresponding one in seq_2 and there are no unknown bases nearby we are able to classify the mutation
      if seq_2[i] != ref[i] and seq_2[i] != 'X' and seq_2[i-2] != 'X' and seq_2[i-1] != 'X' and seq_2[i+1] != 'X' and seq_2[i+2] != 'X':
          
          mut = ""    #mut is the string which describes the mutation according to the nucleotide
          mut_a = ""  #mut_a is the string which describes the mutation according to the aminoacid
          
          #deletion case
          if seq_2[i] == '-':
            mut_type = 4

            mut = mut + 'del ' + str(i+1)
            start = i
            while seq_2[i] == '-' and i < 29673:
              i = i+1
            mut = mut + ':' + str(i)
            end = i

            if (end-start)%3: #if the number of deleted bases isn't a multiple of 3 we have a frame-shift
              mut_type = 6


          #insertion case
          elif ref[i] == '-':
            mut_type = 5
            string_ins = ''
            mut = mut + 'ins ' + str(i+1)
            start = i
            while ref[i] == '-' and i < 29673:
              string_ins = string_ins + seq_2[i]  #evaluating all the inserted bases
              i = i+1
              
            mut = mut + ' ' + string_ins
            end = i

            if (end-start)%3:  #if the number of inserted bases isn't a multiple of 3 we have a frame-shift
              mut_type = 7

          #base substitution case
          else:

            #saving the triplet of bases according to the position in the variable amm
            pos_amm = (i - 265)%3            
            ref_amm = gen_code[ref[i-pos_amm:i+3-pos_amm]]  #ref[i-pos_amm:i+3-pos_amm]   represents the current triplet of reference
            new_triplet = seq_2[i-pos_amm:i+3-pos_amm]      #seq_2[i-pos_amm:i+3-pos_amm] represents the current triplet of seq_2
            amm = gen_code[new_triplet]     

            #silent case
            if ref_amm == amm:
              mut_type = 0

            #nonsense case
            elif ref_amm == 'STOP':
              mut_type = 3
            elif amm == 'STOP':
              mut_type = 1
            
            #missense conservative case
            
            elif amm_class[ref_amm] == amm_class[amm]:
              mut_type = 2

            #missense non-conservative case
            else:
              mut_type = 3
          

          #According to the position and the type of the mutations, the correct counter is updated

          #ORF1ab : 266-21555
          #S (Spike) : 21563-25384
          #ORF3A : 25393-26220
          #E : 26245-26472
          #M : 26523-27191
          #ORF6 : 27202-27387
          #ORF7a : 27394-27759
          #ORF7b : 27760-27887
          #ORF8 : 27894-28259
          #N : 28274-29533
          #ORF10 : 29558-29674

          if i > 264 and i < 21555:
            gene = 'ORF1ab'
            if mut_type == 7:
              ORF1ab[6] = ORF1ab[6] + 1
            else:
              ORF1ab[mut_type] = ORF1ab[mut_type] + 1
            pos = i - 264
          elif i > 21561 and i < 25384:
            gene = 'S'
            if mut_type == 7:
              S[6] = S[6] + 1
            else:
              S[mut_type] = S[mut_type] + 1
            pos = i - 21561
          elif i > 25391 and i < 26220:
            gene = 'ORF3a'
            if mut_type == 7:
              ORF3a[6] = ORF3a[6] + 1
            else:
              ORF3a[mut_type] = ORF3a[mut_type] + 1
            pos = i - 25391
          elif i > 26243 and i < 26472:
            gene = 'E'
            if mut_type == 7:
              E[6] = E[6] + 1
            else:
              E[mut_type] = E[mut_type] + 1
            pos = i - 26243
          elif i > 26521 and i < 27191:
            gene = 'M'
            if mut_type == 7:
              M[6] = M[6] + 1
            else:
              M[mut_type] = M[mut_type] + 1
            pos = i - 26521
          elif i > 27200 and i < 27387:
            gene = 'ORF6'
            if mut_type == 7:
              ORF6[6] = ORF6[6] + 1
            else:
              ORF6[mut_type] = ORF6[mut_type] + 1
            pos = i - 27200
          elif i > 27392 and i < 27759:
            gene = 'ORF7a'
            if mut_type == 7:
              ORF7a[6] = ORF7a[6] + 1
            else:
              ORF7a[mut_type] = ORF7a[mut_type] + 1
            pos = i - 27392
          elif i > 27758 and i < 27887:
            gene = 'ORF7b'
            if mut_type == 7:
              ORF7b[6] = ORF7b[6] + 1
            else:
              ORF7b[mut_type] = ORF7b[mut_type] + 1
            pos = i - 27758
          elif i > 27892 and i < 28259:
            gene = 'ORF8'
            if mut_type == 7:
              ORF8[6] = ORF8[6] + 1
            else:
              ORF8[mut_type] = ORF8[mut_type] + 1
            pos = i - 27892
          elif i > 28272 and i < 29533:
            gene = 'N'
            if mut_type == 7:
              N[6] = N[6] + 1
            else:
              N[mut_type] = N[mut_type] + 1
            pos = i - 28272
          elif i > 29556 and i < 29674:
            gene = 'ORF10'
            if mut_type == 7:
              ORF10[6] = ORF10[6] + 1
            else:
              ORF10[mut_type] = ORF10[mut_type] + 1
            pos = i - 29556
          elif i < 260 or i > 29700:
            if mut_type == 7:
              ORF1ab[6] = ORF1ab[6]
            else:
              ORF1ab[mut_type] = ORF1ab[mut_type]
            pos = -1
          else:
            gene = 'non-coding region'
            pos = -1
            if mut_type == 7:
              NON_COD[6] = NON_COD[6] + 1
            else:
              NON_COD[mut_type] = NON_COD[mut_type] + 1

          #Descriptions of the mutation in the substitution case
          if mut_type < 4:
            mut = mut + ref[i] + str(i+1) + seq_2[i]
            if pos > 0:   #if it's in a coding region
              mut_a = mut_a + ref_amm + str(int(np.ceil(pos/3))) + amm
            else:
              mut_a = '-'
          else: #if it's not substitution
            mut_a = '-'
          
          #row of the mutation list
          r = (mut, mut_a, mutation_types[mut_type], gene)
          mutation_list.append(r)

      #in the end, update the iterator
      i = i+1
  
  #by concatenating all the counters we obtain a row of the statistics dataset
  label = [-1] #label for unknwon variant
  total_row = ORF1ab+S+ORF3a+E+M+ORF6+ORF7a+ORF7b+ORF8+N+ORF10+NON_COD + label
  return mutation_list, total_row

#Function which receives a fasta as input and outputs the dataset of statistics, the dataset with the strings aligned and the list of mutations

def align_and_process(input, string_length = 1000):

  new_data = read(input)                        #sequence to be aligned
  ref=read('original_covid_genome.fasta')[0]    #reference sequence

  #Initialization of the outputs
  col = ["Reference aligned", "Full sequence", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "label"]
  total_aligned = pd.DataFrame(columns = col)
  data=[]
  mutation_list=[]

  
  i=0
  for line in new_data:
    i = i+1
    seq1,seq2,seq_s = align(ref,line, span = string_length, step = int(string_length/2))  #alignment
    m_list,new_row=process_mutations(seq1,seq2)                                                      #processing of the alignment

    data.append(new_row)
    mutation_list=mutation_list+m_list
    total_aligned = total_aligned.append({'Reference aligned': seq1, 'Full sequence': seq2, 'ORF1ab': seq2[265:21555], 'S': seq2[21562:25384], 'ORF3a': seq2[25392:26220], 'E': seq2[26244:26472], 'M' : seq2[26522:27191], 'ORF6' : seq2[27201:27387], 'ORF7a' : seq2[27393:27759], 'ORF7b' : seq2[27755:27887], 'ORF8' : seq2[27893:28259], 'N' : seq2[28273:29533], 'ORF10' : seq2[29557:29674], 'label' : -1}, ignore_index=True)
    
    #verbose
    if not (i+1)%5:
      print('Row ' + str(i+1) + ': ' + str((i+1)/len(new_data)*100) + '%')

  #preparation of the header for the stats dataset
  col_ORF1ab = ["s_ORF1ab", "ns_ORF1ab", "mc_ORF1ab", "mnc_ORF1ab", "del_ORF1ab", "ins_ORF1ab", "fs_ORF1ab"]
  col_S = ["s_S", "ns_S", "mc_S", "mnc_S", "del_S", "ins_S", "fs_S"]
  col_ORF3a = ["s_ORF3a", "ns_ORF3a", "mc_ORF3a", "mnc_ORF3a", "del_ORF3a", "ins_ORF3a", "fs_ORF3a"]
  col_E = ["s_E", "ns_E", "mc_E", "mnc_E", "del_E", "ins_E", "fs_E"]
  col_M = ["s_M", "ns_M", "mc_M", "mnc_M", "del_M", "ins_M", "fs_M"]
  col_ORF6 = ["s_ORF6", "ns_ORF6", "mc_ORF6", "mnc_ORF6", "del_ORF6", "ins_ORF6", "fs_ORF6"]
  col_ORF7a = ["s_ORF7a", "ns_ORF7a", "mc_ORF7a", "mnc_ORF7a", "del_ORF7a", "ins_ORF7a", "fs_ORF7a"]
  col_ORF7b = ["s_ORF7b", "ns_ORF7b", "mc_ORF7b", "mnc_ORF7b", "del_ORF7b", "ins_ORF7b", "fs_ORF7b"]
  col_ORF8b = ["s_ORF8b", "ns_ORF8b", "mc_ORF8b", "mnc_ORF8b", "del_ORF8b", "ins_ORF8b", "fs_ORF8b"]
  col_N = ["s_N", "ns_N", "mc_N", "mnc_N", "del_N", "ins_N", "fs_N"]
  col_ORF10 = ["s_ORF10", "ns_ORF10", "mc_ORF10", "mnc_ORF10", "del_ORF10", "ins_ORF10", "fs_ORF10"]
  col_NON_COD = ["s_NON_COD", "ns_NON_COD", "mc_NON_COD", "mnc_NON_COD", "del_NON_COD", "ins_NON_COD", "fs_NON_COD"]
  label = ["label"]
  col = col_ORF1ab + col_S + col_ORF3a + col_E + col_M + col_ORF6 + col_ORF7a + col_ORF7b + col_ORF8b + col_N + col_ORF10 + col_NON_COD + label
  dataset_stats = pd.DataFrame(data, columns = col) 


  return dataset_stats, total_aligned, mutation_list

#Description of the variant by mutations frequency
def key_mutations(mutation_list, length):
  most_freq_mut = []
  for i,j in enumerate(mutation_list):
      mutation_list[i] = tuple(j)
  #By using listToFreq and sortFreqDict we can sort the mutations by their frequency in the sample
  mutation_freq = listToFreq(mutation_list)
  mutation_freq = sortFreqDict(mutation_freq)
  for s in mutation_freq: 
    if s[0] > length*0.30: 
      
      most_freq_mut.append([s[1][0],s[1][1],s[1][2],s[1][3],s[0]/length*100])  #we append to the output only the ones which occurs at least in 50% of the cases

  key_mutations = pd.DataFrame(most_freq_mut, columns = ['Mutation (Nucleotide)', 'Mutation (Amminoacid)', 'Type', 'Gene','Percentage'])
  return key_mutations

#Function to save the outputs as csv given a name
def export_data(stats, genome_aligned, key_mutations, name = "new"):
  stats.to_csv(name + "_stats.csv", index = False)
  genome_aligned.to_csv(name + "_genomes_aligned.csv", index = False)
  key_mutations.to_csv(name + "_key_mutations.csv", index = False)
