#!/usr/bin/env python
#

from Bio.PDB import *
import numpy as np
import json
import sys
from Bio.PDB import Entity

pdb=sys.argv[1] #input pdb
chain= sys.argv[2] #chain of the residue of interest
lig = sys.argv[3] #residue from which we will get the center of mass
distance= float(sys.argv[4]) #distance to 3D point
distance2= float(sys.argv[5])
parser = PDBParser()
io=PDBIO()
st = parser.get_structure(pdb[:-4],pdb) #get structure from the pdb

rs_dic={} #dictionary used for checking that we use complete residues (with all atoms that compose each one)
for rs in st.get_residues(): #iterate over the residues of the structure
    count=0
    for at in rs.get_atoms(): #iterate over the atoms of each residue
        count+=1 #count how many atoms each residue contain
    rs_dic[rs]=count #dictionary composed of the residue name and the number of atoms each residue contain, e.g. "ASN104":4

res_at=[] #list used to store all atoms of residue of interest
for ch in st.get_chains(): 
    for res in ch.get_residues():
        if res.get_resname() in lig and ch.get_id() in chain:
            for at in res.get_atoms():
                res_at.append(at)

l=[] #list used to store all atoms in structure
for x in st.get_atoms(): #iterate over atoms of structure
    l.append(x)
neighbor=NeighborSearch(l) #set the NeighborSearch function for all atoms in structure

residues=[]
res_to_fix=[]
dyn_res_list=[] #list used to store ALL residues with all atoms at a given distance from the center
fix_res_list=[]
for at in res_at:
    dyn_atoms=neighbor.search(at.get_coord(),distance,level="A") #select all atoms within a given distance from all atoms in residue of interest
    fix_atoms=neighbor.search(at.get_coord(),distance2,level="A") #select all atoms within a given distance from all atoms in residue of interest

    for res in rs_dic: 
        count1=0
        for at in dyn_atoms:
            if at.get_parent()==res: #search for the residues of all the atoms
                count1+=1
        if count1==rs_dic[res]:
            if res not in dyn_res_list:
                dyn_res_list.append(res) #list used to store all residues with ALL the atoms within the distance from the input coordinates

        count2=0
        for at in fix_atoms:
            if at.get_parent()==res: #search for the residues of all the atoms
                count2+=1
        if count2==rs_dic[res]:
            if res not in fix_res_list:
                fix_res_list.append(res) #list used to store all residues with ALL the atoms within the distance from the input coordinates

for res in dyn_res_list: #store all residues encountered from each atom in a residue list, without repeating
    if res in residues:
        residues=residues
    else:
        residues.append(res)
        

for res in fix_res_list: #store all residues encountered from each atom in a residue list, without repeating
    if res in residues:
        residues=residues
    else:
        residues.append(res)
        res_to_fix.append(res)
        
fix=[]
res_to_fix.sort()
for res in res_to_fix:
    if res.get_full_id()[2]==" ":
        fix.append("_"+":"+str(res.get_id()[1]))
    elif res.get_full_id()[2]!="Z":
        fix.append(str(res.get_full_id()[2])+":"+str(res.get_id()[1]))
        
    
print(json.dumps(fix))

class ResSelec(Select): #class used to select the residues and store them in a new pdb file
    def accept_residue(self,residue):
        if residue in residues:
           return True
        else:
            return False

io.set_structure(st)
io.save(pdb[:-4]+"_mod"+str(int(distance))+"_wl_fix.pdb", ResSelec())