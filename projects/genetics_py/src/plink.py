

import logging
import pandas as pd


logger = logging.getLogger(__name__)


def load_fam(fplink_fam):
    cols = 'fid,iid,N1,N2,sex,phe'.split(',')
    fam = pd.read_csv(fplink_fam, sep='\\s+',
                      header=None,
                      dtype={'fid': str, 'iid': str, 'sex': int, 'phe': int},
                      names=cols)
    #fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)
    return fam


def load_fam_plink(fplink, chrom_split=False):
    if chrom_split:
        fplink = fplink.replace('%', '1')
    fplink_fam = fplink + '.fam'
    return load_fam(fplink_fam)


def load_sample_id(fplink, chrom_split=False):
    fam = load_fam_plink(fplink, chrom_split)
    cols = ['fid', 'iid']
    sample_id = fam.loc[:, cols].copy()
    return sample_id
