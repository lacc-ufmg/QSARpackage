import os
import sys
import click
from LQTAQSAR.MD.md import *

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols', '-m',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing molecules'
              )
@click.option('--extension', '-e',
              metavar='<extension>',
              nargs=1,
              required=True,
              help='Extension of the input files'
              )
@click.option('--alignment_file', '-a',
              metavar='<alignment_file>',
              type=click.Path(),
              required=False,
              help='CSV File containing the information about the atoms for alignment'
              )
@click.option('--smarts', '-s',
              metavar='<smarts>',
              nargs=1,
              required=False,
              help='Smarts string to be used in the alignment'
              )
@click.option('--atom', '-p',
              metavar='[atom]',
              multiple=True,
              required=False,
              help='Probe.'
              )
@click.option('--step', '-d',
              metavar='<step>',
              type=float,
              nargs=1,
              required=False,
              help='Steps for navegation on matrix.'
              )

def run(mols,extension,alignment_file,smarts,atom,step):
	md = MolecularDynamics(mols,extension,alignment_file,smarts,atom,step)
	md.run_MD_set()
	if alignment_file is None and smarts is None:
		pass
	else:
		md.runAlignment()
	md.prepare_LQTAGrid_files()
	print(atom)
	print(step)
	if not atom is None and not step is None:
		print("Vai rodar LQTAgrid")
		md.runLQTAGrid()	

if __name__ == '__main__':
	run()