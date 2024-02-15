import argparse

class runlqtagrid:
    """Eu criei isso aqui só pra testar o argparse.
    O ideal é ter uma 'API' que chama as funções de dentro do pacote.
    Aí isso aqui deve ser substituído em breve.
    Mantendo aqui só pra poder fazer o commit."""
    def main(self):
        print(f"\033[1;32mYou've called LQTAGrid, but I'm not implemented yet.\033[0m")
        raise NotImplementedError("This function is not implemented yet.")

class runmd:
    def main(self):
        print(f"\033[1;32mYou've called MD, but I'm not implemented yet.\033[0m")
        raise NotImplementedError("This function is not implemented yet.")


def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=lambda _: parser.print_help())
    subparsers = parser.add_subparsers()

    # LQTAGrid CLI
    parser_lqtagrid = subparsers.add_parser('lqtagrid', help='Command for running LQTAGrid.')
    parser_lqtagrid.set_defaults(func=runlqtagrid.main)
    parser_lqtagrid.add_argument('--mols', '-m', metavar='<path>', required=True, help='Path to the directory containing the molecule files. Both .gro and .itp files should be present.')
    parser_lqtagrid.add_argument('--extension', '-e', metavar='<extension>', required=True, help='Extension of the input molecule files. This should match the file type of the molecules.')
    parser_lqtagrid.add_argument('--coordinates', '-c', metavar='<x> <y> <z>', type=float, nargs=3, help='Coordinates of the box in which the molecules are placed. Specify as three space-separated values.')
    parser_lqtagrid.add_argument('--dimensions', '-d', metavar='<x> <y> <z>', type=float, nargs=3, help='Dimensions of the box in which the molecules are placed. Specify as three space-separated values.')
    parser_lqtagrid.add_argument('--atom', '-a', metavar='[atom]', action='append', required=True, help='Specify the probe atom. This option can be used multiple times to specify multiple atoms.')
    parser_lqtagrid.add_argument('--step', '-s', metavar='<x>', type=float, required=True, help='Step size for navigating the matrix. This determines the granularity of the grid.')
    parser_lqtagrid.add_argument('--output', '-o', metavar='<path_output>', required=True, help='Path to the output file where the matrix will be saved.')

    # MD CLI
    parser_md = subparsers.add_parser('md', help='Command for running MD.')
    parser_md.set_defaults(func=runmd.main)
    parser_md.add_argument('--mols', '-m', metavar='<path>', required=True, type=str, help='Path to the directory containing the molecule files. Both .gro and .itp files should be present.')
    parser_md.add_argument('--extension', '-e', metavar='<extension>', required=True, type=str, help='Extension of the input molecule files. This should match the file type of the molecules.')
    parser_md.add_argument('--alignment_file', '-a', metavar='<alignment_file>', required=False, type=str, help='Path to a CSV file containing information about the atoms for alignment. This is optional.')
    parser_md.add_argument('--smarts', '-s', metavar='<smarts>', required=False, type=str, help='SMARTS string to be used in the alignment. This is optional.')
    parser_md.add_argument('--atom', '-p', metavar='[atom]', action='append', required=False, help='Specify the probe atom. This option can be used multiple times to specify multiple atoms. This is optional.')
    parser_md.add_argument('--step', '-d', metavar='<step>', required=False, type=float, help='Step size for navigating the matrix. This determines the granularity of the grid. This is optional.')
  
  
    args = parser.parse_args()
    args.func(args)