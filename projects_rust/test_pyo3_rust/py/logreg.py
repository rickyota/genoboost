

import numpy as np

from sklearn.linear_model import LogisticRegression


def fib(n):
    a, b = 1, 0
    for _ in range(n):
        a, b = b, a + b
    return b


def fib5():
    return fib(5)


def ar_sum(ar):
    return ar.sum()


def ar_twice(ar):
    # return np.array([0.0, 1.0])
    return ar * 2


# covs: N x M
def logreg(covs, y, sample_weight):
    # print(covs.shape)
    # print(y.shape)
    # print(sample_weight.shape)

    # fit_intercept=False

    clf = LogisticRegression(penalty='none', random_state=0).fit(covs, y, sample_weight)
    intercept = clf.intercept_
    coef = clf.coef_
    # print(intercept,coef)
    return np.concatenate([intercept, coef[0]], 0)
