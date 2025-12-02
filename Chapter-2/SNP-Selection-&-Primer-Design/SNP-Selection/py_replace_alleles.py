#A.Jacquemart 
#20 March 2023
#!/usr/bin/ env python
import sys
import re

#To run, use the command python3 py_replace_alleles.py name of fasta file.fa name of allele info file.txt
fasta = open(sys.argv[1],'r').readlines() #open OG fasta file
alleles= open(sys.argv[2],'r').readlines() #open the allele info file 
output=open("name_of_output.fa",'w') #create output file

#splitting apart the fasta file into relevant components
headers=[] #creates list of all of the >NCXXXX ids
seqs=[]#creates a list of all of genomic sequences
for line in fasta: #goes through each line of the fasta file
	if re.search(">NC_",line):
		line=line.replace("\n","")
		line=line.replace(">","")
		headers.append(line) #adds >NC line to list if that's in the line
	else:
		line=line.replace("\n","")
		seqs.append(line) #otherwise, saves the genomic sequences to the other list

#create a dictionary to associate each ">NC" id to it's fasta sequence
pairs={} 
for header in headers: #loops through the list of ">NC" ids
	index=headers.index(header)
	pairs[header]=seqs[index] #adds the ">NC" id as the key and the genomic sequence as the value
#print(pairs)

#create another dictionary for associate the ">NC" ids to their allele combos 
SNPs={}
for line in alleles: #loops through each line of the file
	line=line.split("\t") #splits the lines by tabs
	SNPs[line[0]]=line[1].replace("\n","") #adds the ">NC" ids as key and the allele combo as the value
#print(SNPs)		

#then, generates the output file
for header in headers: #loops through the ">NC" ids
	var_seq= pairs[header] #saves the associated genomic sequence to a new variable
	bp=len(var_seq) # determines how long the fasta sequence is and saves to new variable
	front=var_seq[0:100] #get substring of just the first set of characters up until the SNP of interest in the genomic sequence. NOTE: python indexes start at 0 
	back=var_seq[101:] #get substring of the last characters of the sequence after the N representing the SNP of interest
	SNP_alleles=SNPs[header] #get the allele combo associated to the ">NC" id
	output.write(">"+str(header)+"\n"+ str(front)+ str(SNP_alleles)+str(back)+"\n") #write these pieces to a new fasta file
output.close() #close output file	
