# -*- coding: utf-8 -*-
"""


@author: Shirly
"""
'''
The script finds all the amino acids in the interface between protein A and protein B in the complex.
This is a complementary script to the script: ComplexContactPipe.py
'''
cutoff = 4.0

from Bio.PDB import *
from Bio.PDB import PDBParser


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

def get_all_atoms(residue_map):
	all_atoms_out = []
	for residue in residue_map:
		for atom in residue_map[residue]:
			all_atoms_out.append(atom)
			atom.set_bfactor(0.0)
	return all_atoms_out


def is_contact(res_1,other_atoms,cutoff):
	
	 for atom in res_1:
		  ns = NeighborSearch(other_atoms)
		  center = atom.get_coord()
		  neighbors = ns.search(center, cutoff) # 5.0 for distance in angstrom
		  residue_list = Selection.unfold_entities(neighbors, 'R') # R for residues
		  if len(residue_list)>0:
			  return True
	 return False 

def get_contacts(struc,all_atoms,verbose,cutoff):
    progress = 0 
    contacts = {}
    for residue in struc:
        progress += 1
        
        
           
        atom_list = struc[residue]
        #print("atom_list:",atom_list)
        outcome = is_contact(atom_list,all_atoms,cutoff)
        if outcome:
            if (not residue.split("_")[0] in contacts.keys()):
                contacts[residue.split("_")[0]] = []
            contacts[residue.split("_")[0]].append(residue)
    return contacts





def find_residues(pdbID,first_chain,second_chain):
    pdbFile = 'complexes/' + pdbID + '.pdb'
    str_1 = PDBParser().get_structure('first_one', pdbFile)
    str_2 = PDBParser().get_structure('second_one', pdbFile)
    atoms_1 = Selection.unfold_entities(str_1, 'C')
    atoms_2 = Selection.unfold_entities(str_2, 'C')
    input_1 = get_atom_list(atoms_1,first_chain)
    input_2 = get_atom_list(atoms_2,second_chain)
    all_atoms_1 = get_all_atoms(input_1)
    all_atoms_2 = get_all_atoms(input_2)
    contacts_1 = get_contacts(input_1,all_atoms_2,"First molecule, residue ",cutoff)
    contacts_2 = get_contacts(input_2,all_atoms_1,"Second molecule, residue ",cutoff)
    count=0
    count1=0
    count2=0
    for c in first_chain:
        if c in contacts_1.keys():
            for residue in contacts_1[c]:
                count+=1
                count1+=1
    for c in second_chain:
        if c in contacts_2.keys():
            for residue in contacts_2[c]:
                count+=1
                count2+=1
    return count1,count2,count
        


if __name__ == '__main__':
    data=open("data.txt","r")
    results=open("aa_in_contact_4.0.txt","w")
    results.write("pdb_ID"+"\t"+"aa in first chain"+"\t"+"aa in second chain"+"\t"+"total\n")
    for line in data:
        line1=line.split()
        pdbID=line1[0][:4]
        print(pdbID)
        first_chain=line1[0].split("_")[1].split(":")[0]
        print(first_chain)
        second_chain=line1[0].split("_")[1].split(":")[1]
        print(second_chain)
        first,second,num_of_residues=find_residues(pdbID,first_chain,second_chain)
        print(num_of_residues)
        pdbAndChains=pdbID+"_"+first_chain+":"+second_chain+"\t"
        #line=line[:-1]
        #results.write(line+"\t"+str(num_of_residues)+"\n")
        results.write(pdbAndChains+"\t"+str(first)+"\t"+str(second)+"\t"+str(num_of_residues)+"\n")
    results.close() 
    data.close()
    print("Done\n")