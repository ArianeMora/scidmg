###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import argparse
import os
import sys

from sciutil import SciUtil

from scidmg import __version__
from scidmg import SciDMG


def print_help():
    lines = ['-h Print help information.']
    print('\n'.join(lines))


def run(args):
    dmg = SciDMG(args.dmr, args.pm, args.dmc,
               dmc_methdiff=args.c_d, dmc_padj=args.c_p, dmc_uid=args.c_id, dmc_padj_cutoff=args.c_c,
               dmr_methdiff=args.r_d, dmr_padj=args.r_p, dmr_uid=args.r_id, dmr_padj_cutoff=args.r_c,
               dmg_filename=args.o, min_perc_agreement=args.mp, min_meth_diff=args.md, plot=args.v
            )
    dmg.run()
    dmg.print_stats()


def gen_parser():
    parser = argparse.ArgumentParser(description='scie2g')
    parser.add_argument('--dmr', type=str, help='DMR file output from a program like MethylSig (csv) annotated to genes.')
    parser.add_argument('--dmc', type=str, help='DMC file output from a program like MethylKit (csv) annotated to genes.')
    parser.add_argument('--g', type=str, help='Gene ID column (must exist in all files passed). Rename the columns '
                                              'to contain the same ID if it is not the same column name.')
    # Optional.
    parser.add_argument('--pm', type=str, help='Percentage of methylation, annotated to genes (optional).')
    parser.add_argument('--o', type=str, default='DMG_output.csv', help='Output file (full path)')
    parser.add_argument('--c_id', type=str, default='uid', help='Unique identifier for the CpG in the DMC file.')
    parser.add_argument('--c_d', type=str, default='meth.diff', help='Column in DMC file with the change in DNA methylation.')
    parser.add_argument('--c_p', type=str, default='qvalue', help='Column in DMC file with the p adjusted value.')
    parser.add_argument('--c_c', type=int, default=0.1, help='Minimum padjusted value for the DMCs.')

    parser.add_argument('--r_id', type=str, default='uid', help='Unique identifier for the CpG in the DMR file.')
    parser.add_argument('--r_d', type=str, default='meth_diff', help='Column in DMR file with the change in DNA methylation.')
    parser.add_argument('--r_p', type=str, default='fdr', help='Column in DMR file with the p adjusted value.')
    parser.add_argument('--r_c', type=float, default=0.1, help='Minimum padjusted value for the DMRs.')

    parser.add_argument('--mp', type=float, default=0.6, help='Minimum percentage of agreement between DMCs and DMRs (default is 60% hence 0.6)')
    parser.add_argument('--md', type=float, default=0.01, help='Minimum amount of change in DNA methylation for DMRs (default is 0.01)')

    parser.add_argument('--v', type=bool, default=False, help='Make volcano plots.')

    return parser


def main(args=None):
    parser = gen_parser()
    u = SciUtil()
    if args:
        sys.argv = args
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'scidmg v{__version__}')
        sys.exit(0)
    else:
        print(f'scidmg v{__version__}')
        args = parser.parse_args(args)
        # Validate the input arguments.
        if not os.path.isfile(args.dmr):
            u.err_p([f'The DMR file could not be located, file passed: {args.dmr}'])
            sys.exit(1)
        if not os.path.isfile(args.dmc):
            u.err_p([f'The DMC file could not be located, file passed: {args.dmc}'])
            sys.exit(1)
        # Otherwise we have need successful so we can run the program
        u.dp(['Running sci-DMR on DMR file: ', args.dmr,
              '\nWith DMC file: ', args.dmc,
              '\nPercentage methylation file: ', args.pm,
              '\nSaving to output file: ', args.o,
              '\nGene ID column that exists in both DMC file and DMR file: ', args.g,
              '\nDMC unique ID column:', args.c_id,
              '\nDMC meth diff column:', args.c_d,
              '\nDMC padj column:', args.c_p,
              '\nDMR unique ID column:', args.r_id,
              '\nDMR meth diff column:', args.r_d,
              '\nDMR padj column:', args.r_p,
              '\nMinimum percentage of agreement:', args.mp,
              '\nMinimum methylation difference: ', args.md,
              '\nDMR padj cutoff:', args.r_c,
              '\nDMC padj cutoff: ', args.c_c
              ])
        # RUN!
        run(args)
    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()
    # ----------- Example below -----------------------
    # root_dir = '../'
    # main(["--dmr", root_dir + "tests/data/methylSig_prom.csv", "--dmc", root_dir + "tests/data/methylKit_DMC.csv"])