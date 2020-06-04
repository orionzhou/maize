
def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = ''
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("", help = '',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = x)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()
