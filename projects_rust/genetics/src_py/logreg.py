

import numpy as np
from sklearn.linear_model import LogisticRegression


# covs: N x M
def logreg_sklearn(covs, y, sample_weight):
    # print(covs.shape)
    # print(y.shape)
    # print(sample_weight.shape)

    # fit_intercept=False

    # no randomness
    # no use-> fixed by random_state
    #np.random.seed(51)

    clf = LogisticRegression(penalty=None, random_state=0).fit(covs, y, sample_weight)

    #clf = LogisticRegression(penalty=None, random_state=0, max_iter=1000).fit(covs, y, sample_weight)
    #clf = LogisticRegression(penalty='none', random_state=0, max_iter=1000).fit(covs, y, sample_weight)
    # stay the same
    #clf = LogisticRegression(penalty='none', random_state=0, max_iter=10000).fit(covs, y, sample_weight)
    # max_iter=100 is not enough; different from smartcore
    #clf = LogisticRegression(penalty='none',random_state=0).fit(covs, y,sample_weight)
    #clf = LogisticRegression(penalty='none',random_state=0,solver='newton-cg').fit(covs, y,sample_weight)

    intercept = clf.intercept_
    coef = clf.coef_
    # print(intercept,coef)
    return np.concatenate([intercept, coef[0]], 0)
