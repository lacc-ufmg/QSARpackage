#!/usr/bin/env python
import os
import pandas as pd
import click
from qsar_package.md.align import *

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols', '-m',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing CEPs from MD simulation'
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

def run(mols,alignment_file,smarts):
	if alignment_file is None and smarts is None:
		print("An alignment file or a smarts string should be provided")
	elif not alignment_file is None:
		runAF(mols,alignment_file)
	elif not smarts is None:
		runSMARTS(mols,smarts)

if __name__ == '__main__':
	run()
