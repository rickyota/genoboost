# TODO: use path

import argparse
import logging
import numpy as np
from sklearn.model_selection import KFold, ShuffleSplit

from . import plink
from . import logger_setting

logger = logging.getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cross_validation',
                        action='store_true',
                        help='file to output')

    parser.add_argument('--cross_validation_n',
                        type=int,
                        help='file to output')

    parser.add_argument('--dout',
                        help='dir to output')

    parser.add_argument('--fplink',
                        help='')

    parser.add_argument('--funrelated',
                        help='file of unrelated samples')

    parser.add_argument('--test_prop',
                        type=float,
                        default=0.2,
                        help='prop of test')

    parser.add_argument('--seed',
                        type=int,
                        default=51,
                        help='random seed')

    args = parser.parse_args()
    logger.debug(args)
    return args


def split_obs_ts(sample_id, test_prop, seed=None):
    rs = ShuffleSplit(n_splits=1, test_size=test_prop, random_state=seed)
    (obs_index, test_index) = list(rs.split(sample_id))[0]

    obs_index = np.sort(obs_index)
    obs_id = sample_id.loc[obs_index, :]
    test_index = np.sort(test_index)
    test_id = sample_id.loc[test_index, :]

    return obs_id, test_id


def split_tr_va_cv(sample_id, cv_n, seed=None):
    kf = KFold(n_splits=cv_n, shuffle=True, random_state=seed)

    # [[tr0, va0], ...]
    idxs = list(kf.split(sample_id))
    idxs = [[np.sort(v) for v in vs] for vs in idxs]

    ids = [[sample_id.loc[v, :] for v in vs] for vs in idxs]

    assert len(ids) == cv_n

    return ids


# split into obs and test, then split obs into tr and va in n-fold
# TODO: same case prop
def split_cv(dout, fplink, cv_n, test_prop, seed):
    sample_id = plink.load_sample_id(fplink)

    obs_id, test_id = split_obs_ts(sample_id, test_prop, seed)

    # np.random.seed(seed)
    #whole_n = len(sample_id)
    #test_n = int(whole_n * prop_test)
    #whole_index = sample_id.index.to_numpy()
    #whole_index_perm = np.random.permutation(whole_index)
    #test_index = whole_index_perm[:test_n]
    # test_index.sort()
    #test_id = sample_id.loc[test_index, :]
    #obs_index = whole_index_perm[test_n:]
    # obs_index.sort()
    #obs_id = sample_id.loc[obs_index, :]

    ftest = dout + 'test.samples'
    test_id.to_csv(ftest, sep='\t', index=False)

    fobs = dout + 'obs.samples'
    obs_id.to_csv(fobs, sep='\t', index=False)

    tr_va_ids = split_tr_va_cv(sample_id, cv_n, seed)
    for (cvi, (tr, va)) in enumerate(tr_va_ids):
        ftr = dout + 'tr.cv' + str(cvi) + '.samples'
        tr.to_csv(ftr, sep='\t', index=False)

        fva = dout + 'va.cv' + str(cvi) + '.samples'
        va.to_csv(fva, sep='\t', index=False)


def main():
    logger.debug("main()")

    args = argument()

    if args.cross_validation:
        if args.funrelated is None:
            split_cv(args.dout, args.fplink, args.cross_validation_n, args.test_prop, args.seed)
        else:
            print('extract unrelated in training dataset')
            NotImplementedError()

    logger.debug("Done!!")


if __name__ == '__main__':
    logger_setting.setting()
    main()
