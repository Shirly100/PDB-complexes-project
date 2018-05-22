# -*- coding: utf-8 -*-
"""


@author: malki,Shirly
"""
'''
The script finds in a complex which chains in protein A contact chains in protein B.
The script takes into account only different chains and reduce simillar chains. 
'''
import os
import Bio

import csv
import os
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
import Bio.PDB
from pypdb.pypdb import *
import pypdb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#consts
cutoff = 3.9
i_cutoff = 10.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
def SelectColumn(fileName, column, delim  = "\t"):
    
    """The function select some specific columns from the file"""
    
    #open the file and separate the specific column using a known delim
    with open (fileName, encoding = "utf8") as inf:
        reader = csv.reader (inf, delimiter = delim)
        ReqCol = list(zip(*reader))[column]
        
    return ReqCol


def find_chains_in_protein(pdbEntry):
    '''find_chains_in_protein -> find the chanis are belong to each protein
        return {'entry': [[[chains in protein 1 divided by entities]], [[chains in protein 2 divided by entities]]]
        for example:
            {'1A0O': [[['A', 'C', 'E', 'G']], [['B', 'D', 'F', 'H']]]}
        }'''
    #get the protein's details from th pdb
    print(pdbEntry)
    tmp = pypdb.get_all_info(pdbEntry)['polymer']
    prevProtein = {}
    prevIndex = -1
    protein = {pdbEntry : []}
    #iterate over entity and capsulate into same proteins
    for entity in tmp:
        if (type(entity['chain']) == type(dict())):
            chains = [i for i in entity['chain'].values()]
        else:
            chains = [i['@id'] for i in entity['chain']]
        #find the different attribute between the 2 dictionaries via xor
        unmatched_item = set(entity["polymerDescription"].items()) ^ set(prevProtein.items())
        
        #if the dictionary are indentical so the chains belong to the same protein:
        # if the polymerDescription is not indentical they are different proteins
        if len(unmatched_item) > 0 :
            prevIndex = prevIndex + 1
            prevProtein = entity["polymerDescription"]
            #put the group chains in the aproppreate protein
            protein[pdbEntry].append([])
            protein[pdbEntry][prevIndex].append(chains)
        
        # match to the same protein
        elif len(unmatched_item) == 0 :
            protein[pdbEntry][prevIndex].append(chains)
    
    return protein

def find_duplicates( dic1,file ): #dic1 is in the format {"1EAY":[[['A','B'],['C','D']]]}
    #try:
    
        for i in dic1:#e.g 1EAY 
            pdb_id=len(dic1[i])
            print("pdb_id",pdb_id)
            count=0
            while(count<pdb_id):
                macro=len(dic1[i][count])#number of macroMolecules in the complex (e.g 1EAY->2)
                print(dic1[i][count])#e.g 1EAY->[['A', 'B'], ['C', 'D']]
                dic2={}
                list=[]
                count1=0
                chains=len(dic1[i][count][0])#number of chains in macromolecule e.g 1EAY->2
                print("chains",chains)
                if( count==0 and chains==1):
                    break
                flag=0
                
                for j in range(0, chains):
                    c=(dic1[i][0][0][j])
                    print("c",c)
                    if(c==''):
                        break
                    p = PDBParser()
                    dicD={}
                    
                    
                    structure = p.get_structure('x', file)
                    for model in structure:
                        chain = model[c]
                        for residue in chain:
                            if not("het=W" in ascii(residue)):
                                if not("het=H" in ascii(residue)):
                                    start = ascii(residue).find("resseq") 
                                    end = ascii(residue).find(" ", start)
                                    resseq=ascii(residue)[start:end]
                                    start1=resseq.find("=") 
                                    positionListD=resseq[start1+1:]
                                    dicD[positionListD]=residue.get_resname()
                    answer=True
                    for key in dicD:
                        if(key in dic2):
                            if(dicD[key]!=dic2[key]):
                                    answer=False
                            else:
                                dic2[key]=dicD[key]
                                
                    if(answer==True and j!=0):
                        dic1[i][0][0][j]=""
                    count+=1
            print(dic1)
            return dic1
                #break
        
        
        

def get_atom_list(structure, chains):
    output = dict()
    for chain in structure:
        if chain.id in chains:
            for residue in chain.get_residues():
                hetflag, resseq, icode=residue.get_id()
                the_id = (chain.id+"_"+str(resseq)+"_"+icode).strip()
                for atom in residue.get_unpacked_list():
                    if hetflag==' ':
                        if the_id in output:
                            output[the_id].append(atom)
                        else:
                            output[the_id] = [atom]
    return output   
                 
#checks the contact between residues
def is_contact(res_1,other_atoms,cutoff):
	 #check for every atom in residue the posibble neigbors
	 for atom in res_1:
		  ns = NeighborSearch(other_atoms)
		  center = atom.get_coord()
		  neighbors = ns.search(center, cutoff) # 5.0 for distance in angstrom
		  residue_list = Selection.unfold_entities(neighbors, 'R') # R for residues
		  if len(residue_list)>0:
			  return True
	 return False 	

#return posibble contacts per chain (from getInterface)
def get_contacts(struc,all_atoms,verbose,cutoff):
    progress = 0 
    contacts = {}
    for residue in struc:
        progress += 1
        #if len(verbose)>0:
            #print (verbose,progress,"out of",len(struc))
        atom_list = struc[residue]
        outcome = is_contact(atom_list,all_atoms,cutoff)
        if outcome:
            if (not residue.split("_")[0] in contacts.keys()):
                contacts[residue.split("_")[0]] = []
            contacts[residue.split("_")[0]].append(residue)
    return contacts

             			
#Filter out all the atoms from the chain,residue map given by residue_map
def get_all_atoms(residue_map):
	all_atoms_out = []
	for residue in residue_map:
		for atom in residue_map[residue]:
			all_atoms_out.append(atom)
			#Set the b-factor to zero for coloring by contacts
			atom.set_bfactor(0.0)
	return all_atoms_out

	#return the chains are in contact with the other protein in complex	
def in_contact_chains(str_1, str_2, chains_1, chains_2):
    #Load the structures - they can be the same!
    atoms_1 = Selection.unfold_entities(str_1, 'C') # C for chains
    atoms_2 = Selection.unfold_entities(str_2, 'C') # C for chains
    #get the mapping from chain,residue id to the atom lists
    input_1 = get_atom_list(atoms_1,chains_1)
    input_2 = get_atom_list(atoms_2,chains_2)
    #get the full atom lists for neighbor search
    all_atoms_1 = get_all_atoms(input_1)
    all_atoms_2 = get_all_atoms(input_2)
    contacts_1 = get_contacts(input_1,all_atoms_2,"First molecule, residue ",cutoff)
    contacts_2 = get_contacts(input_2,all_atoms_1,"Second molecule, residue ",cutoff)
    return contacts_1.keys(), contacts_2.keys()
    #run neighbor search on both instances - not optimal but good enough for most imaginable applications.
	 	 

def get_complexes_contacts(fileName):
    #creat a list of complexe's pdb entries
    #pdbEntries = [i.replace("_","") for i in SelectColumn(fileName, 0)]
    #open the complexe's info. file
    inputFile = open(fileName, 'r') 
    #output file always like: with_chains_<file name>
    outputFile = open('contact_chains3_' + fileName, 'w+')
    # for every complex do the next pipeline
    #for complexID in pdbEntries:
        #read the information about the complex
    for line in inputFile:
        
       
    #line = inputFile.readline()
        #find chains per protein in complex
        line1=line.split("\t")
        complexID=line1[0].replace('\t',"")
       
        print(complexID," * ")
        if(len(complexID) is 4):
            chainsPerProtein = find_chains_in_protein(complexID)
            print(chainsPerProtein)
            #get the path to the approppriate complex pdb file
            pdbFile = 'complexes/' + complexID.lower() + '.pdb'
            #reduce indent chains in protein
            reduceDuplicate = find_duplicates(chainsPerProtein, pdbFile)
            #load the complex
            str_1 = PDBParser().get_structure('first_one', pdbFile) # load your molecule
            str_2 = PDBParser().get_structure('second_one', pdbFile) # load your molecule
            #get the cains for each protein
            chains_1 = ''.join([''.join(chain) for chain in reduceDuplicate[complexID][0]])
            chains_2 = ''.join([''.join(chain) for chain in reduceDuplicate[complexID][1]])
            #get the chains in contact per protein in complex
            contact_chains1, contact_chains2 = in_contact_chains(str_1, str_2, chains_1, chains_2)
            # change the complex entry to entry_<protein 1 contact chains>:<protein 2 contact chains>
            tmp = line.split('\t')
            tmp[0] = tmp[0] + '_' + ''.join(contact_chains1) + ':' +''.join(contact_chains2)
            outputFile.write('\t'.join(tmp))
            print(tmp[0])
            #print(reduceDuplicate)
            #print(chains_1, reduceDuplicate[complexID][0], contact_chains1)
            #print(chains_2, reduceDuplicate[complexID][1], contact_chains2)
        
    inputFile.close();
    outputFile.close();
if __name__ == '__main__':
    os.chdir(os.path.dirname(__file__) or os.getcwd())
    #invoke pipeline function
    get_complexes_contacts("data_without_chains.txt")
    print("Done\n")

    
    