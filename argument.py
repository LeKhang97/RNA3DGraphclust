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
					choices = ['G', 'M', 'L', 'C', 'H'],
					help="Clustering algorithm. Either: G (Girvan-Newman); M (Markov); L (Louvain); C (Clauset-Newman-Moore); H (Hierachical-based) (default).")
    
    # Subparser for -a G
    #parser_a_G = subparsers.add_parser('G', help='Arguments for Girvan-Newman clustering algorithm')
    
    # Subparser for -a M
    parser_a_M = subparsers.add_parser('M', help='Arguments for Markov clustering algorithm')
    parser_a_M.add_argument('-e', type= int, default= 5, help='Cluster expansion factor (default = 5)')
    parser_a_M.add_argument('-i', type= float, default= 1.4, help='kernel type (default = 1.4)') 
    
    # Subparser for -a L
    parser_a_L = subparsers.add_parser('L', help='Arguments for Louvain clustering algorithm')
    parser_a_L.add_argument('-r', type=float, default= 1, help='resolution of the algorithm (default = 1)')

    # Subparser for -a C
    parser_a_C = subparsers.add_parser('C', help='Arguments for Clauset-Newman-Moore clustering algorithm')
    parser_a_C.add_argument('-r', type=float, default= 1, help='resolution of the algorithm (default = 1)')

    # Subparser for -a H
    parser_a_H = subparsers.add_parser('H', help='Arguments for Hierarchical clustering algorithm')
    parser_a_H.add_argument('-t', type=float, default= 0.4, help='upper threshold of modularity score for top-down algo (default = 0.4)')
    parser_a_H.add_argument('-b', type=float, default= 0.1, help='ratio of nodes for bottom-up algo to calculate threshold (default = 0.1)')
    parser_a_H.add_argument('-r', type=float, default= 1, help='resolution of the algorithm (default = 1)')
    
    args = parser.parse_args()      
    
    return args

def process_args():
    args = main_argument()
    largs = [args.input, args.algorithm, args.distance, args.weight]
    largs2 = [args.outfile, args.verbose, args.chain]

    algo_list = ['Girvan-Newman', 'Markov', 'Louvain', 'CNM', 'Hierachical-based']
    
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
        if not hasattr(args, 'r'):
            args.r = 1

        if args.verbose:
            print(f"r: {args.r}")
        largs += [args.r]
    
    elif args.algorithm == 'C':
        if not hasattr(args, 'r'):
            args.r = 1

        if args.verbose:
            print(f"r: {args.r}")
        largs += [args.r]
        
    elif args.algorithm  == 'H':
        if not hasattr(args, 't'):
            args.t = 0.4
        if not hasattr(args, 'b'):
            args.b = 0.1
        if not hasattr(args, 'r'):
            args.r = 1
        
        if args.verbose:
            print(f"t: {args.t}, b: {args.b}, r: {args.r}")
        largs += [args.t, args.b, args.r]
            
    else:
        sys.exit("Unrecognized algorithm!")

    largs += [args.threshold]

    return largs, largs2