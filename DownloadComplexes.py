# -*- coding: utf-8 -*-
"""


@author: malki, Shirly
"""
'''
This script is running on linux environment.
The script downloads complexes proteins from the PDB sites and creates new PDB files, containing only 
the chains in contact between protein A and protein B of the complex.
'''
from time import sleep
import os
import csv
import urllib
def SelectColumn(fileName, column, delim  = "\t"):
    
    """The function select some specific columns from the file"""
    
    #open the file and separate the specific column using a known delim
    with open (fileName) as inf:
        reader = csv.reader (inf, delimiter = delim)
        ReqCol = list(zip(*reader))[column]
    return ReqCol

def download_complex_file(complexes):
    #the maximum amount of files to zip file 
    index = 5
    length = len (complexes) 
    counter = 0
    message = "\"" + str(index) + " tgz file was created \""
    #tmp directory
    os.system('rm -Rf complexes')
    os.system('mkdir complexes')
    
    os.system('mkdir complexes_chains_in_contact')
    #download the complexes pdb file in the file
    for c in complexes:
         #download and save the pdb file in tmp
	 urllib.urlretrieve('https://files.rcsb.org/download/' + c.replace('_','') + '.pdb', 'complexes/' + c.lower() + '.pdb')
	 print "Downloading " + c + "..."
         #pdbl.retrieve_pdb_file(c, pdir='tmp')
         
         
     
         counter = counter + 1
         if counter == length:
             #compress every $length files
             os.system("tar -czvf complexes_all_chains" + ".tgz complexes")
       
             length = 0
             sleep(30)
    data=open('data_with_contact_chains_3.9.txt',"r")
    for line in data:
        list_chains=[]
        line=line.split('\t')
        pdbID=line[0][0:4]
        print(pdbID)
        chains_in_contact=line[0].split("_")[1]
        print(chains_in_contact)
        chains=chains_in_contact.split(":")
        list_chains.append( chains[0])
        seconds=chains[1]
        for i in range(0, len(seconds)):
            c=seconds[i]
            list_chains.append(c)
        print(list_chains)
        pdbFile = 'complexes/' + pdbID.lower() + '.pdb'
        #os.system('cd complexes_chains_in contact')
        #newPdbFile =  pdbID.lower() + '.pdb'
        newPdbFile='complexes_chains_in_contact/' +pdbID.lower() + '.pdb'
        new=open(newPdbFile ,"w")
        filepdb=open(pdbFile ,"r")
        for line in filepdb:
            line1 = line.split()
            if not(line1[0] == 'ATOM'):
                new.write(line)
            elif line1[0] == 'ATOM':
                cha= line1[4]
                #print(cha)
                if  cha in list_chains :
                    new.write(line)
        new.close

    

if __name__ == '__main__':
    os.chdir(os.path.dirname(__file__) or os.getcwd())
    #the complexes files
    fileName = "data_without_chains.txt"
    #pdb entry column
    column = 0
    #the delim separates the file
    delim = "\t"
    try:
    	ComplexList = SelectColumn(fileName, column, delim  = "\t")
    	download_complex_file(ComplexList)
    except StandardError:
    	print("Done\n")

    
    