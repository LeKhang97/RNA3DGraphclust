import argparse
from pathlib import Path
from Functions import *

def main_argument():
    parser = argparse.ArgumentParser(description ="This tool is used to detect RNA domains from given 3D coordinate :))")

    subparsers = parser.add_subparsers(dest='subcommand')

    parser.add_argument('-v', '--verbose',
        action ='store_true', 
        help ='verbose mode.')
    
    parser.add_argument('-i',
        '--input',
        required=True,
        help ='input file. Must be in pdb format.')
    
    parser.add_argument('-c',
        '--chain',
        action='store_true',
        help ='process all chains at once. If not, the program will process each chain individually.')
    
    parser.add_argument('-t',
        '--threshold',
        default= 30,
        type= int,
        help ='Lower threshold for sequence length')

    parser.add_argument('-o', '--outfile',
                #default = None,
                action ='store',
                help ='output file.')
    
    parser.add_argument('-d', '--distance',
                        default = 8,
                        type = float,
                        help = 'distance threshold to determine contact between two residues (default = 8)')
    
    parser.add_argument('-w', '--weight',
                        default = "False",
                        type = str,
                        help = 'weight the contact matrix (default = False). If True, the contact matrix will be weighted using the sigmoid function.')

    parser.add_argument('-a', 
					'--algorithm',
                    default = 'G',
					choices = ['G', 'M', 'L', 'H'],
					help="Clustering algorithm. Either: G (Girvan-Newman); M (Markov); L (Louvain); H (Hierachical-based) (default).")
    
    # Subparser for -a G
    #parser_a_G = subparsers.add_parser('G', help='Arguments for Girvan-Newman clustering algorithm')
    
    # Subparser for -a M
    parser_a_M = subparsers.add_parser('M', help='Arguments for Markov clustering algorithm')
    parser_a_M.add_argument('-e', type= int, default= 5, help='Cluster expansion factor (default = 5)')
    parser_a_M.add_argument('-i', type= float, default= 1.4, help='kernel type (default = 1.4)') 
    
    # Subparser for -a L
    parser_a_L = subparsers.add_parser('L', help='Arguments for Louvain clustering algorithm')
    parser_a_L.add_argument('-n', type=int, default= 2, help='number of clusters (default = 2)')

    # Subparser for -a H
    parser_a_H = subparsers.add_parser('H', help='Arguments for Hierarchical clustering algorithm')
    parser_a_H.add_argument('-m', type=float, default= 0.4, help='upper threshold of modularity score for top-down algo (default = 0.4)')

    args = parser.parse_args()      
    
    return args

def process_args():
    args = main_argument()
    largs = [args.input, args.algorithm, args.distance, args.weight]
    largs2 = [args.outfile, args.verbose, args.chain]

    algo_list = ['Girvan-Newman', 'Markov', 'Louvain', 'Hierachical-based']
    
    algo = [i for i in algo_list if i[0] == args.algorithm][0]

    if args.verbose:
        if args.weight in ("True", "T", "true"):
            print("Using weighted matrix")
        else:
            print("Using non-weighted matrix")

        print("Using algorithm: ", algo)
        print(f"Arguments for {algo}:")
        
    if args.algorithm == 'G':
        if args.verbose:
            print(f"This algorithm does not require any arguments.")

    elif args.algorithm == 'M':
        if not hasattr(args, 'e'):
            args.e = 5
        
        if not hasattr(args, 'i'):
            args.i = 1.4

        if args.verbose:
            print(f"e: {args.e}, i: {args.i}")
        largs += [args.e, args.i]

    elif args.algorithm == 'L':
        if not hasattr(args, 'n'):
            args.n = 2

        if args.verbose:
            print(f"n: {args.n}")
        largs += [args.n]
        
    elif args.algorithm  == 'H':
        if not hasattr(args, 'm'):
            args.m = 0.4
        
        if args.verbose:
            print(f"m: {args.m}")
        largs += [args.m]
            
    else:
        sys.exit("Unrecognized algorithm!")

    largs += [args.threshold]

    return largs, largs2
