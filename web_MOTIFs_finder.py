#!/usr/bin/python3

"""
Title: web_MOTIF_finder
    
Author: Maria Tsalkitzidou

Description:
    This script uses Flask in order to build a web application. This application scans a set of fasta sequences (submitted by the user) for a provided motif or multiple motifs (submitted by the user as a file or manually in a textbox). 
    It produces a tab-separated table with the sequence ids of the sequences that have the motif(s), the number of times the motif(s) exist in each sequence and its position in the sequence. The application also produces a new fasta file (one fasta file per motif) with only the sequences that have the motif(s).

List of functions:
    1. check_fastafile()
    2. get_fasta_in_list()
    3. get_list_from_file()
    4. get_motifs()
    5. accepted_nuc()

List of modules:
    1. flask
    2. shutil
    3. re
    4. os
    
List of bugs:
    1. The lines of the table are not aligned properly
    2. The code is slow with large datasets
    3. The web application crushes with very big datasets
    4. If the user provides a fasta file with mixed DNA and protein sequences the program will not notify him of the problem and it will run.
    
     
    
Usage:
     python web_MOTIFs_finder.py
     Ctrl+click to open the web page.
     Ctrl+c to stop the program
"""

#%% Import packages

import re
import os
import shutil
from flask import Flask, render_template, request

app = Flask(__name__)



#%% Functions -----------------------------------

# Function for checking if the fasta file has headers and if every header has a sequence and vice versa
def check_fastafile(lid, lseq):
    #If headers were found    
    if len(lid) >= 1:        
        #Check if the number of headers is different than the number of sequences.
        if len(lid) != len(lseq):
            message='The amount of IDs and sequences should be the same. Please make sure the input file is a fasta file'                    
        else:
            message='OK'                
    # If no headers were found
    else:
        message='The sequences must have their ID as header. Please make sure the input file is a fasta file'
        
    return(message)
        
#Function to get fasta sequences in lists
def get_fasta_in_list(fasta_file):    
    lid=[]
    lseq=[]
    temp=''
    
    #Iterate through the fasta file
    for line in fasta_file:
        line=line.decode('utf-8') # change the code of upload file
        line=line.strip()
        
        if line.startswith(' '):
            message='The fasta file is wrong.'
            return render_template('error_handling.html', message=message)
        
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
        #if the file is not empty
        accepted_nuc(temp) #check that the sequence is DNA, RNA or protein
        lseq.append(temp)
        
        
    #Check if the fasta file is in the appropriate format.
    if check_fastafile(lid, lseq)=='OK':
        return(lseq, lid)
    else:
        message= check_fastafile(lid, lseq)
        return render_template('error_handling.html', message=message)

#Function for checking if the sequences contain only DNA, RNA, protein and N
def accepted_nuc(sequence):
    
    #Define two lists with the accepted letters
    acc_nuc = ["A","T","G","C", "U", "N"]
    acc_prot = ["F", "L", "S", "Y", "C", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "G", "D", "E", "*"]
    
    #Iterate though the sequence
    for nuc in sequence:
        #Check if the sequence has DNA, RNA or protein
        if nuc not in acc_nuc or nuc not in acc_prot:
            message= 'The sequence appears to not be DNA, RNA or protein.'
            return render_template('error_handling.html', message=message)
        else:
            continue
        
    return    

#Function to store data from column in file as a list
def get_list_from_file(filename):
    mylist=[]
    for line in filename:
        line=line.decode('utf-8')
        line=line.upper()
        l=line.strip()
        mylist.append(l)
        
    return(mylist)


#Function to get the motif occurences
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
            seq_dict[sequence_id]=sequence
            
    
    return(motif_dict, seq_dict)

#%% Preprocessing of file format -----------------------------

# limitation of upload file format.
fasta_format = ['fasta','fna','fa', 'faa']


#%% Web ------------------------------------------

#Define the submission page
@app.route('/')
def motif_finder_input():
    return render_template('motif_finder_input.html')

#Define the results page
@app.route('/', methods=['POST'])
def show_motif_finder():
    
    # specifiy the HTTP method to use when sending form-data
    if request.method == 'POST':
    
        #For manually entering the motif(s)
        motif=request.form['motif']
         
        
        #When the user entered a file with motifs
        motif_file=request.files['motif_file']
        
        #Access the fasta file
        fasta_file=request.files['fastafile']
               
        #Check if the user provided a fasta file
        if fasta_file.filename == '': #if no file was provided
            message='No fasta file selected'
            return render_template('error_handling.html', message=message)                   
        
        #If there is already a FASTA_files directory from a previous run, delete it and create a new empty one to store the results frm the current run
        if os.listdir('./static/FASTA_files'):
            shutil.rmtree('./static/FASTA_files')
            os.makedirs('./static/FASTA_files')
        
                                

        if fasta_file: #check the file and file format
            #Read the fasta file and store the sequences and their ids in lists, with the function get_fasta_in_list
            fasta_seq=get_fasta_in_list(fasta_file)
            #Assign the results of the function in two variables
            seq_ids=fasta_seq[1]
            seq=fasta_seq[0]
                            
            
            #Get the number of scanned sequences
            number_seq=len(seq)
            
            #Create a header for the results table
            Header=['#SequenceID', 'Abundance', 'Position']
            
            # Empty variables to store the results
            motif_output = {}
            motif_name =[]
            ab_motifs={}
            
            
        #
            #Open an output file
            output = open('./static/motif_info.tsv','w')
            
            #If the user provided motifs manually
            if motif:
                motif_list=re.split(';|\n', motif) #Split the multiple motif(s)
                if motif_file:# prevent conflicts if users submit sequence and file at the same time
                    message='Enter the motif(s) manually in the textbox OR in a file. Not both.'
                    return render_template('error_handling.html', message=message)                   
                
                else:
                    #Iterate through the entered motifs to find matches for each motif
                    for mo in motif_list:                                                       
                        mo=mo.upper().strip() #Make the motifs in uppercase to avoid confilcts
                        
                        #Use the function to get the motifs and assign the results in new variables 
                        results=get_motifs(mo, seq, seq_ids)
                        motif_results=results[0]
                        sequence_results=results[1]
        #
                        #Add the motif to the output list of motif names
                        motif_name.append(mo)
                        
                        # Add the first key in the dictionary
                        motif_output[mo] = []
                        
                        #Create a counter for the number of sequences having the motif
                        counter=0
                        
                        #Add the header to the output file
                        output.write('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))
                        
                        #Print the results in an appropriate way
                        for k in motif_results: #Iterate through the results
                            if len(motif_results[k]) >=1: #If there was at least one motif match in the sequence
                                counter+=1 #There was a match so add 1 to the counter
                                for p in motif_results[k]: #Add the sequence ID, the abundance of the motif in that sequence and its position to the results dictionary                                           
                                    motif_output[mo].append([k, len(motif_results[k]), p]) #add the reults to the dictionary
                                    output.write('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p)) #print output to the file
                        #Add the number of sequences having at least one motif match to the results dictionary
                        ab_motifs[mo]=counter
                        
                        
                        #Open a file named by the motif 
                        f=open('./static/FASTA_files/{}.fasta'.format(mo), 'w')                                       
                        
                        #Print the sequences that have the motif in a new fasta file        
                        for l in sequence_results:           
                            f.write('{}\n{}\n'.format(l, sequence_results[l]))
                                                           
                        
                        #Close the file
                        f.close()
                        
        #
            
            #If the user provided a file of one or multiple motifs
            elif motif_file:
                if motif:# prevent conflicts if users submit sequence and file at the same time
                    message='Enter the motif(s) manually in the textbox OR in a file. Not both.'
                    return render_template('error_handling.html', message=message)                   
                else:
                    #Use the function to get the motifs from the file
                    motif=get_list_from_file(motif_file)                   
                    #Check that the file was not empty
                    if len(motif) >= 1:
                        #If it is not empty iterate through the motifs to find matches for each motif
                        for mo in motif:
                            #Use the function to get the motifs and assign the results in new variables
                            results = get_motifs(mo, seq, seq_ids)
                            motif_results = results[0]
                            sequence_results=results[1]
                            
                            # Add the motif to the output list of motif names
                            motif_name.append(mo)
                            
                            # Add the first key in the dictionary
                            motif_output[mo] = []                           
                            
                            #Create a counter for the number of sequences having the motif
                            counter = 0                           
        #
                            #Add the header to the output file
                            output.write('##Motif:{}\n#SequenceID\tAbundance\tPosition\n'.format(mo))                           
                            
                            # Print output to website & file
                            for k in motif_results: #Iterate through the results
                                if len(motif_results[k]) >= 1: #If there was at least one motif match in the sequence
                                    counter += 1 #There was a match so add 1 to the counter
                                    for p in motif_results[k]: #Add the sequence ID, the abundance of the motif in that sequence and its position to the results dictionary
                                        motif_output[mo].append([k, len(motif_results[k]), p]) #add the reults to the dictionary
                                        output.write('{}\t{}\t{}\n'.format(k, len(motif_results[k]), p)) #print output to the file
                            #Add the number of sequences having at least one motif match to the results dictionary
                            ab_motifs[mo]=counter
                            
                            
                    
                            #Open a file named by the motif 
                            f=open('./static/FASTA_files/{}.fasta'.format(mo), 'w')                                       
                            
                            #Print the sequences that have the motif in a new fasta file        
                            for l in sequence_results:           
                                f.write('{}\n{}\n'.format(l, sequence_results[l]))
                                                               
        #
                            #Close the file
                            f.close()
                            
                    else: #If the user gave a blank motif file lead him to the error page
                        message='The motif file is empty'
                        return render_template('error_handling.html', message=message)
            
            #Close the motif info file
            output.close()
            
            #Store all the fasta files in a zip folder
            shutil.make_archive('./static/FASTA_files', 'zip', './static/FASTA_files')
                    
                
            return render_template('motif_finder.html', number_seq=number_seq, ab_motifs=ab_motifs, motif_name=motif_name, header=Header, motif_output=motif_output)
        


        
        else: #If the fasta file is empty lead the user to the error page
            message='The fasta file is empty'
            return render_template('error_handling.html', message=message)


if __name__ == '__main__':
    app.debug = True
    app.run()