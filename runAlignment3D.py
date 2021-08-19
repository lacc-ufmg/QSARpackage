from openbabel import pybel
import os
import sys
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
from rdkit.Chem.Draw.IPythonConsole import drawMol3D
from openbabel import openbabel
import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols', '-m',
              metavar='<mols>',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing the molecules to align'
              )
@click.option('--extension','-e',
              metavar='<extension>',
              required=True,
              help='Extension of yhe input files'
              )
@click.option('--alignment_file', '-a',
              metavar='<alignment_file>',
              type=click.Path(),
              required=True,
              help='CSV File containing the information about the atoms for alignment'
              )
@click.option('--output_dir', '-o',
              metavar='<output_dir>',
              nargs=1,
              required=True,
              help='Directory where aligned molecules will be saved'
              )


def run(mols,alignment_file,extension,output_dir):

    df = pd.read_csv(os.path.join(mols,alignment_file),sep=";",header = None)

    format_file = "g09" if (extension == "log" or extension == "out") else extension

    obc = openbabel.OBConversion()
    obc.SetOutFormat("mol2")

    # Pybel
    mols = [list(pybel.readfile(format_file,os.path.join(mols,f)))[0]
            for f in df[0]]

    # refMol = Chem.MolFromMol2File(os.path.join(input_dir,df[0][0]),removeHs=False)

    refMol = Chem.MolFromMolBlock(mols[0].write("mol"),removeHs=False)

    atoms = [[int(a)-1 for a in df[1][j].split(',')] for j in range(len(df[1]))]

    refAtoms = atoms[0]

    m = pybel.readstring("mol",Chem.MolToMolBlock(refMol))
    charges = [a.partialcharge for a in mols[0]]
    for i,a in enumerate(openbabel.OBMolAtomIter(m.OBMol)):
        a.GetPartialCharge()
        a.SetPartialCharge(charges[i])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filename = df[0][0].split('.')[0]+".mol2"
    obc.WriteFile(m.OBMol,os.path.join(output_dir,filename))

    for i in range(1,len(mols)):
        print(df[0][i])
        # prbMol = Chem.MolFromMol2File(os.path.join(input_dir,df[0][i]),removeHs=False)
        prbMol = Chem.MolFromMolBlock(mols[i].write("mol"),removeHs=False)
        prbAtoms = atoms[i]
        alinhamento = AlignMol(prbMol,refMol,atomMap=list(zip(prbAtoms,refAtoms)))
        m = pybel.readstring("mol",Chem.MolToMolBlock(prbMol))
        charges = [a.partialcharge for a in mols[i]]
        for j,a in enumerate(openbabel.OBMolAtomIter(m.OBMol)):
            a.GetPartialCharge()
            a.SetPartialCharge(charges[j])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        filename = df[0][i].split('.')[0]+".mol2"            
        obc.WriteFile(m.OBMol,os.path.join(output_dir,filename))

if __name__ == '__main__':
    # input_dir = sys.argv[1]
    # alignment_file = sys.argv[2]
    # output_dir = sys.argv[3]
    # run(input_dir,alignment_file,output_dir)
    run()
