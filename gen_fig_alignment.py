from openbabel import pybel
import os
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions #Only needed if modifying defaults
from rdkit.Chem.Draw.IPythonConsole import drawMol3D
from rdkit.Chem import rdFMCS

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def gen_fig_alignment(directory:str, extension:str, format_file:str, output:str, smartsString:str)->None:

    # directory = "/home/jpam/google_drive/UFMG/luiz_claudio/tatiane/moleculas/gjf/mol2"
    # extension = "mol2"
    # format_file = "mol2"

    # Pybel
    mols = [list(pybel.readfile(format_file,os.path.join(directory,f)))[0] 
            for f in os.listdir(directory) if f.endswith(extension)]

    str_mols = ""
    for m in mols:
        str_mols += m.write('smi')
        str_mols += '\n'

    # Rdkit
    rdmols = []
    for m in mols:
        rdmols.append(Chem.MolFromPDBBlock(m.write('pdb'),removeHs=False))
    
    if smartsString == "" or smartsString == None:
        res=rdFMCS.FindMCS(rdmols)
        smartsString = res.smartsString

    # Prepare molecules labeling atoms with numbers and generating 2D depictions
    for m in rdmols:
        mol_with_atom_index(m)
        AllChem.Compute2DCoords(m)

    common_atoms = [m.GetSubstructMatch(Chem.MolFromSmarts(smartsString)) for m in rdmols]

    DrawingOptions.includeAtomNumbers=True
    img = Draw.MolsToGridImage(rdmols, 
                            molsPerRow=4, 
                            useSVG=True,subImgSize=(500, 500), 
                            legends=[m.title for m in mols],
                            highlightAtomLists = common_atoms,
                            maxMols=110,
                            )

    with open(os.path.join(directory,output),'w') as f:
        f.write(img.data)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Generate figure alignment')
    parser.add_argument('-d','--directory', type=str, help='Directory with molecules')
    parser.add_argument('-e','--extension', type=str, help='Extension of the files',default='mol2')
    parser.add_argument('-ff','--format_file', type=str, help='Format of the files',default='mol2')
    parser.add_argument('-o','--output', type=str, help='Output file',default='fig.svg')
    parser.add_argument('-ss','--smartsString', type=str, help='SMARTS string',required=False)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    gen_fig_alignment(args.directory, args.extension, args.format_file, args.output, args.smartsString)