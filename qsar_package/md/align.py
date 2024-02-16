from openbabel import pybel
import os
import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolAlign import AlignMol

def alinhaPAC(mol,atoms,filename):

    refMol = mol[0]

    for m in mol[1:]:
        alinhamento = AlignMol(m,refMol,atomMap=list(zip(atoms,atoms)))

def alinhaPACs(refMol,refAtoms,prbMol,prbAtoms,filename):

    w = Chem.SDWriter(filename)

    for mol in prbMol:
        alinhamento = AlignMol(mol,refMol,atomMap=list(zip(prbAtoms,refAtoms)))
        w.write(mol)
    w.close()

def alinhaPAC_SMARTS(mol,patt):
    refMol = mol[0]
    atoms = list(refMol.GetSubstructMatch(patt))
    for m in mol[1:]:
        alinhamento = AlignMol(m,refMol,atomMap=list(zip(atoms,atoms)))

def alinhaPACs_SMARTS(refMol,prbMol,patt,filename):

    refAtoms = list(refMol.GetSubstructMatch(patt))
    prbAtoms = list(prbMol[0].GetSubstructMatch(patt))

    w = Chem.SDWriter(filename)

    for mol in prbMol:
        alinhamento = AlignMol(mol,refMol,atomMap=list(zip(prbAtoms,refAtoms)))
        w.write(mol)
    w.close()

def convertPDB(filename):
	os.system('obabel {} {} {}'.format(
	    filename,
	    '-O',
	    filename[:-3]+"pdb"
	))
	os.remove(filename)

def runAF(directory,alignment_file):

	df = pd.read_csv(os.path.join(directory,alignment_file),sep=";",header = None)

	files = [f for f in os.listdir(directory) if f.endswith("pdb")]

	files.sort()

	mols = [pybel.readfile('pdb',os.path.join(directory,f))
	    for f in files if f.endswith("pdb")]

	rdmols = []
	for mol in mols:
		print(mol)
		rdconf = []
		for c in mol:
			print(c)
			rdconf.append(Chem.MolFromPDBBlock(c.write('pdb'),removeHs=False))
		rdmols.append(rdconf)

	refAtoms = [int(a)-1 for a in df[1][0].split(',')]
	refmol = rdmols[0][0]


	for i,mol in enumerate(rdmols):
		parts = re.search("(.*)\.(\w+)",files[i])
		name = parts.group(1)
		atoms = [int(a)-1 for a in df[1][i].split(',')]
		filename = os.path.join(directory,name+"_aligned.mol")
		alinhaPAC(mol,atoms,filename)
		alinhaPACs(refmol,refAtoms,mol,atoms,filename)
		convertPDB(filename)

def runSMARTS(directory,smarts):

	patt = Chem.MolFromSmarts(smarts)

	files = [f for f in os.listdir(directory) if f.endswith("pdb")]

	files.sort()

	mols = [pybel.readfile('pdb',os.path.join(directory,f))
	    for f in files if f.endswith("pdb")]

	rdmols = []
	for mol in mols:
	    rdconf = []
	    for c in mol:
	        rdconf.append(Chem.MolFromPDBBlock(c.write('pdb'),removeHs=False))
	    rdmols.append(rdconf)

	refmol = rdmols[0][0]

	for i,mol in enumerate(rdmols):
		parts = re.search("(.*)\.(\w+)",files[i])
		name = parts.group(1)
		filename = os.path.join(directory,name+"_aligned.mol")
		alinhaPAC_SMARTS(mol,patt)
		alinhaPACs_SMARTS(refmol,mol,patt,filename)
		convertPDB(filename)
