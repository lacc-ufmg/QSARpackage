#!/usr/bin/env
# coding: utf-8

import click
from qsar_package.lqtagrid.grid_generate import *
#from . import grid_generate


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols','-m',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='files path, gro and itp.'
              )
@click.option('--extension','-e',
              metavar='<extension>',
              required=True,
              help='Extension of yhe input files'
              )
@click.option('--coordinates', '-c',
              metavar='<x> <y> <z>',
              type=float,
              nargs=3,
              required=False,
              help='Coordinates of the box.'
              )
@click.option('--dimensions', '-d',
              metavar='<x> <y> <z>',
              type=float,
              nargs=3,
              required=False,
              help='Dimensions of the box.'
              )
@click.option('--atom', '-a',
              metavar='[atom]',
              multiple=True,
              required=True,
              help='Probe.'
              )
@click.option('--step', '-s',
              metavar='<x>',
              type=float,
              nargs=1,
              required=True,
              help='Steps for navegation on matrix.'
              )
@click.option('--output', '-o',
              metavar='<path_output>',
              type=click.Path(),
              required=True,
              help='Output matrix file.'
              )
def run(mols, extension, coordinates, dimensions, atom, step, output):
    '''LQTAgridPy is a python version of LQTAgrid, a practical application of
    4D QSAR analysis methodology developed at University of Campinas.

    More: https://github.com/rougeth/LQTAgridPy
    '''
    grid = GridGenerate(
        extension,
        coordinates,
        dimensions,
        atom,
        mols,
        step
    )
    grid.saveGrid(output)


if __name__ == '__main__':
    run()
