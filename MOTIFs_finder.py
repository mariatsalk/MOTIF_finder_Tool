#!/usr/bin/python3
"""
Title: MOTIF_finder
    
Author: Maria Tsalkitzidou

Description:
    This script takes as input a fasta file and a motif or multiple motifs and searches for this motif(s) in  the fasta file. It informs the user about the sequences the motif(s) was/were found in, the number of times the motif(s) exist in each sequence and its position in the sequence. The results are printed on the terminal but they can also be outputted in a tab-separated file. The program can also extract the sequences that have the motif(s) in a new fasta file, for further analyisis.

List of functions:
    1. check_fastafile()
    2. get_fasta_in_list()
    3. get_list_from_file()
    4. get_motifs()
    5. accepted_nuc()

List of modules:
    1. argparse
    2. path form pathlib
    3. re
    4. os
    
List of bugs:
    1. The columns of the table are not equally distanced
    2. The code is slow with large datasets
    3. If the user provides a fasta file with mixed DNA and protein sequences the program will not notify him of the problem and it will run.
    
    
Usage: python MOTIFs_finder --file example_realfasta.fasta --motif GGTGGT --motif_file realmotif.txt --outSeq --outMotif motif_table.tsv
    
 
"""


#%% Import Packages

import argparse
from pathlib import Path
import re
import os

#%% Importing the file from the terminal
parser = argparse.ArgumentParser()

usage ='Find motif(s) in your sequences. Provide a fasta file with sequences and your desired motif or motifs. The motif(s) can be provided manually or in a file using the corresponding flag. The results contain the sequence ids of the sequences that have the motif alongside with the number of times the motifs was found in the sequence and its position. The results are printed and they can also be saved in a file if the appropriate flag is used. The sequences that have the motif can be extracted in a new fasta file if the appropriate flag is used.'

parser.add_argument(
    '--file',
    dest='fasta_file',
    required=True,
    help='Fasta file containing the sequences to search the motif(s).')
parser.add_argument(
    '--motif',
    required=False,
    nargs='+',
    help='Provide one or multiple motifs manually for discovery. For multiple motifs separate them with a semicolon (;).')
parser.add_argument(
    '--motif_file',
    required=False,
    help='Provide one of multiple motifs in a file for discovery. For mutliple motifs they have to be in a column')
parser.add_argument( 
    '--outSeq',
    required=False,
    dest='out_file1',
    action='store_true',
    help='Fasta file with ONLY the sequences that have the motif(s).')
parser.add_argument( 
    '--outMotif',
    dest='out_file',
    required=False,  
    help='A tab-separated file with the motif information.')

args=parser.parse_args()

#%% Functions

# Function for checking if the fasta file has headers and if every header has a sequence and vice versa
def check_fastafile(lid, lseq):
    #If headers were found    
    if len(lid) != 0:
        #Check if the number of headers is different than the number of sequences.
        if len(lid) != len(lseq):
            raise Exception("The amount of IDs and sequences should be the same. Please make sure the input file is a fasta file")
        
        else:
            return True
    
    # If no headers were found
    else:
        raise Exception("The sequences must have their ID as header. Please make sure the input file is a fasta file")

# Function to read through a fasta file and store the sequence ids and the correspronding sequences in two lists.
def get_fasta_in_list(fasta_file):    
    lid=[]
    lseq=[]
    temp=''
     
    for line in fasta_file:
        # Search for the header
        if line.startswith('>'):
            IDline=line.rstrip() #To remove any characters at the end a string
            ID=IDline.split()[0] # To keep the id from the header         
            lid.append(ID) # Adding the ids in the lid list                           
            #print(ID) #Checkpoint            
            
            if temp:
                accepted_nuc(temp) #check that the sequence is DNA, RNA or protein
                lseq.append(temp); #Keep the previous sequence in seq
            temp='' #zero temp
        
        else: # If it found sequence            
            temp += line.rstrip().upper() #remove new line character and make the sequence uppercase because it might be mixed lowercase uppercase
            
    if temp:
        accepted_nuc(temp) #check that the sequence is DNA, RNA or protein        
        #if the file is not empty
        lseq.append(temp)
    
    
    #Check if the fasta file is in the appropriate format.
    if check_fastafile(lid, lseq):
        return(lseq, lid)
    else:
        raise Exception("The fasta file is not in the appropriate format.")
        
        
#Function for checking if the sequences contain only DNA, RNA, protein and N
def accepted_nuc(sequence):
    
    #Define two lists with the accepted letters
    acc_nuc = ["A","T","G","C", "U", "N"]
    acc_prot = ["F", "L", "S", "Y", "C", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "G", "D", "E", "*"]
    
    #Iterate though the sequence
    for nuc in sequence:
        #Check if the sequence has DNA, RNA or protein
        if nuc not in acc_nuc or nuc not in acc_prot:
            raise Exception('The sequence appears to not be DNA, RNA or protein.')
            
        else:
            continue
        
    return

#Function to store data from column in file as a list
def get_list_from_file(filename):
    mylist=[]
    for line in filename:
        l=line.rstrip()
        mylist.append(l)
        
    return(mylist)


#Function to get the motif occurences and the sequences that contain them
def get_motifs(motif, seq, seq_ids):
        
    #Create empty dictionaries to store the results
    motif_dict={}
    seq_dict={}
    
    #Create an regex object with the motif
    pat=re.compile(motif)
    
    #Find the motif instances in the sequences
    for s in seq:
        m_id=seq_ids[seq.index(s)]
        motif_dict[m_id]=[]
        for m in pat.finditer(s):
            motif_dict[m_id].append(m.span())               
       
    #Extract the sequences that have the motif
    for x in range(len(seq)):
        sequence=seq[x]
        if motif in sequence:
            sequence_id=seq_ids[seq.index(sequence)]
            #seq_dict[sequence_id]=''
            seq_dict[sequence_id]=sequence
    
    return(motif_dict, seq_dict)



#%% Load the files

motif=[]
seq=[]
seq_ids=[]

#Check if all the mandatory information information is passed:
if args.motif or args.motif_file and args.fasta_file:
    
    #Open the fasta file
    infile=open(args.fasta_file, 'r')
    
    # Checking if the fasta file is empty
    if os.stat(args.fasta_file).st_size == 0:
        raise Exception("The fasta file should not be empty")
    
    #Store the fasta sequences and their ids in lists using the function
    fasta_seq=get_fasta_in_list(infile)
    seq_ids=fasta_seq[1]
    seq=fasta_seq[0]
    
    #Empty list to store the motif names
    motif_name=[]
    
    #check if motif is passed manually or in a file
    #If it's passed manually assign it to a variable
    if args.motif:
        motif=args.motif
        if args.motif_file: #If the user also provided a motif file
            raise Exception("Provide the motif(s) manually OR in a file. Not both.")
        else:
            print('----------Job processing----------') #Inform the user about the status of the job
            
            #Iterate through every motif
            for mo in motif:                                                       
                mo=mo.upper().strip() #change the motif letters to upper
                
                #Use the function to get the motifs and assign the results in new variables 
                results=get_motifs(mo, seq, seq_ids)
                motif_results=results[0]
                sequence_results=results[1]
                
                #Add the motif to the output list of motif names
                motif_name.append(mo)
                
                #Print the results
                #Print the header
                print('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))
                #Print motif information
                for k in motif_results:
                    if len(motif_results[k]) >=1: #If there was a match in the search
                        for p in motif_results[k]: #Print the table                                       
                            print('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p))
                
                #If the user wants to store the table in a file
                if args.out_file:
                    #Open the output file
                    motif_output=open('{}.tsv'.format(args.out_file), 'a')
                    
                    #Add the header
                    motif_output.write('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))
                    #Print output to file
                    for k in motif_results:
                        if len(motif_results[k]) >=1: #If there was a match in the search
                            for p in motif_results[k]:                                           
                                motif_output.write('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p))
                                
                    #Close the file
                    motif_output.close()
                
                #If the user wants to extract the sequences that have the motif from the original dataset
                if args.out_file1:                   
                    #Open a file named by the motif 
                    seq_output=open('{}.fasta'.format(mo), 'w')
                                
                    #Print the sequences that have the motif in a new fasta file        
                    for l in sequence_results:           
                        seq_output.write('{}\n{}\n'.format(l, sequence_results[l]))
                        
                    #Close the file
                    seq_output.close()
                    
                    
            
        
    #If it's passed in a file check if the file exists and then assign it to a variable
    elif args.motif_file:
        #Check if the file exists
        if Path(args.motif_file).is_file():
            motif_file=open(args.motif_file, 'r')
            motif=get_list_from_file(motif_file)
        else:
            raise Exception("The file does not exist.")
            
        if args.motif: #If the user also entered a motif manually
            raise Exception("Provide the motif(s) manually OR in a file. Not both.")
        else:
            print('----------Job processing----------') #Inform the user about the status of the job
            for mo in motif:                                                       
                mo=mo.upper().strip() #change the motif letters to upper
                
                #Use the function to get the motifs and assign the results in new variables
                results=get_motifs(mo, seq, seq_ids)
                motif_results=results[0]
                sequence_results=results[1]
                
                #Add the motif to the output list of motif names
                motif_name.append(mo)
                
                #Print the results
                #Print the header
                print('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))
                #Print motif information
                for k in motif_results:
                    if len(motif_results[k]) >=1: #If there was a match in the search
                        for p in motif_results[k]:                                       
                            print('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p))
                
                #If the user wants an output file with the motif information
                if args.out_file:
                    motif_output=open('{}.tsv'.format(args.out_file), 'a')
                    
                    #Add the header
                    motif_output.write('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))
                    #Print output to file
                    for k in motif_results:
                        if len(motif_results[k]) >=1: #If there was a match in the search
                            for p in motif_results[k]:                                           
                                motif_output.write('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p))
                                
                    #Close the file
                    motif_output.close()
                
                #If the user wants to extract the sequences that have the motif
                if args.out_file1:
                    #Open a file named by the motif 
                    seq_output=open('{}.fasta'.format(mo), 'w')
                                
                    #Print the sequences that have the motif in a new fasta file        
                    for l in sequence_results:           
                        seq_output.write('{}\n{}\n'.format(l, sequence_results[l]))
                        
                    #Close the file
                    seq_output.close()
                    
    print("----------Job finished----------") #Inform the user about the status of the job
    
    #Close the fasta file
    infile.close()
    
else: #If the user didn't provide the mandatory arguments
    raise Exception("Provide the correct arguments.")