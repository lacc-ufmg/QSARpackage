#!/usr/bin/env
# coding: utf-8

import click
from LQTAQSAR.LQTAGrid import hull_generate
#from . import grid_generate


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--mols','-m',
              metavar='<path>',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing mol files.'
              )
@click.option('--extension','-e',
              metavar='<extension>',
              required=True,
              help='Extension of yhe input files'
              )
@click.option('--atom', '-a',
              metavar='[atom]',
              multiple=True,
              required=True,
              help='Probe.'
              )
@click.option('--initial', '-i',
              metavar='<initial>',
              type=float,
              nargs=1,
              required=False,
              help='Initial distance from molecule hull.'
              )
@click.option('--angle', '-d',
              metavar='<angle>',
              type=float,
              nargs=1,
              required=False,
              help='Angle step.'
              )
@click.option('--radius', '-r',
              metavar='<radius>',
              type=float,
              nargs=1,
              required=False,
              help='Radius step.'
              )
@click.option('--layers', '-l',
              metavar='<layers>',
              type=int,
              nargs=1,
              required=False,
              help='Number of layers.'
              )
@click.option('--output', '-o',
              metavar='<path_output>',
              type=click.Path(),
              required=True,
              help='Output matrix file.'
              )
def run(mols, extension, atom, angle, initial, radius, layers, output):
    '''LQTAgridPy is a python version of LQTAgrid, a practical application of
    4D QSAR analysis methodology developed at University of Campinas.

    More: https://github.com/rougeth/LQTAgridPy
    '''
    hull = hull_generate.HullGenerate(extension, atom, mols, angle, initial, radius, layers)
    hull.saveGrid(output)


if __name__ == '__main__':
    run()
