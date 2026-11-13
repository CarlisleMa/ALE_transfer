import pickle
from openmm.app import PDBFile
with open('/task/output/protein_final_state.pkl','rb') as f: top,pos=pickle.load(f)
with open('/task/output/protein_for_packmol.pdb','w') as f:
 PDBFile.writeFile(top,pos,f,keepIds=True)
print('WROTE',len(list(top.residues())),len(list(top.atoms())))
