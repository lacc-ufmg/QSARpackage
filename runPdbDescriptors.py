import os
from openbabel import pybel
import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--directory','-d',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='files path, gro and itp.'
              )
@click.option('--matrix', '-m',
              metavar='<matrix>',
              type=click.Path(),
              required=True,
              help='CSV File containing the matrix with descriptors'
              )              
@click.option('--input_file', '-i',
              metavar='<input_file>',
              type=click.Path(),
              required=True,
              help='Input molecule file'
              )              
@click.option('--output_file', '-o',
              metavar='<output_file>',
              type=click.Path(),
              required=True,
              help='Output molecule file'
              )


def geraPdbDescritores(directory,matrix,input_file,output_file):
    f = open(os.path.join(directory,matrix),"r")
    descritores = f.readline()
    d = descritores.split(';')
    coord = []
    for v in d[1:]:
        coord.append(v.split('_'))
    coordenadas = []
    for c in coord:
        coordenadas.append((float(c[0]),float(c[1]),float(c[2])))
    f.close()

    extension = input_file.split('.')[-1]
    format_file = "g09" if (extension == "log" or extension == "out") else extension

    if format_file != "pdb":
        mol = list(pybel.readfile(format_file,os.path.join(directory,input_file)))[0]
        pdb = mol.write("pdb")
    else:
        f = open(os.path.join(directory,input_file),'r')
        pdb = f.read()
        f.close()
    f = open(os.path.join(directory,output_file),'w')
    f.write(pdb)
    f.write("TITLE     Descritor\n")
    f.write("MODEL        2\n")
    for i in range(len(coordenadas)):
        f.write("ATOM{:6d}  X   DES     1      {:6.3f}  {:6.3f}  {:6.3f}  1.00  0.00\n".format(i+1,coordenadas[i][0],coordenadas[i][1],coordenadas[i][2]))
    f.write("TER\n")
    f.write("ENDMDL\n")
    f.close()

if __name__=='__main__':
    geraPdbDescritores()
