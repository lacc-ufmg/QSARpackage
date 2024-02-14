import argparse
# from . import run, runmd

class runlqtagrid:
    def main(self):
        print("lqtagrid")

class runmd:
    def main(self):
        print("md")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    parser_lqtagrid = subparsers.add_parser('lqtagrid')
    parser_lqtagrid.set_defaults(func=runlqtagrid.main)

    parser_md = subparsers.add_parser('md')
    parser_md.set_defaults(func=runmd.main)

    args = parser.parse_args()
    args.func(args)