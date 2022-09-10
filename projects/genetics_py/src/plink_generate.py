

import argparse
import logging
import numpy as np

from . import plink
from . import logger_setting

logger = logging.getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cov',
                        action='store_true',
                        help="generate cov")

    parser.add_argument('--pheno',
                        action='store_true',
                        help="generate phe")

    parser.add_argument('--fout',
                        default=None,
                        help='file to output, if None, fplink')

    parser.add_argument('--fplink',
                        help='prefix of plink file')

    parser.add_argument('--seed',
                        type=int,
                        default=51,
                        help='random seed')

    args = parser.parse_args()
    logger.debug(args)
    return args


# assume no chrom split
# sex: 0: male, 1: female
def generate_cov(fout, fplink, seed):
    np.random.seed(seed)

    fam = plink.load_fam_plink(fplink)
    cov = fam[['fid', 'iid']].copy()

    if '0' in fam['sex']:
        raise RuntimeError('Unknown sex.')
    sex = 2 - fam['sex'].astype(int)
    cov.loc[:, 'sex'] = sex

    age = np.random.randint(50, 70, size=len(cov))
    cov.loc[:, 'age'] = age

    cols = ['fid', 'iid', 'sex', 'age']
    cov = cov.loc[:, cols]

    fcov = fout + '.cov'
    #fcov = fplink + '.cov'
    cov.to_csv(fcov, sep='\t', index=False)


# assume no chrom split
# phe: 1: control, 2: case
def generate_phe(fout, fplink, seed):
    np.random.seed(seed)

    fam = plink.load_fam_plink(fplink)
    phe = fam[['fid', 'iid']].copy()

    case_prop = 0.3
    phe_bi = np.random.choice(2, size=len(phe), p=[1 - case_prop, case_prop])
    phe_bi = 1 + phe_bi
    phe.loc[:, 'phe'] = phe_bi

    cols = ['fid', 'iid', 'phe']
    phe = phe.loc[:, cols]

    fcov = fout + '.pheno'
    #fcov = fplink + '.cov'
    phe.to_csv(fcov, sep='\t', index=False, header=False)


def main():
    logger.debug("main()")

    args = argument()

    if args.fout is not None:
        fout = args.fout
    else:
        fout = args.fplink

    if args.cov:
        generate_cov(fout, args.fplink, args.seed)

    if args.pheno:
        generate_phe(fout, args.fplink, args.seed)

    logger.debug("Done!!")


if __name__ == '__main__':
    logger_setting.setting()

    main()
