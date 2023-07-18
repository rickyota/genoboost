

#import numpy as np
#import statsmodels.api as sm


""" def calc_likelihood_regression_score(score, phe):
    # llmodel= score + const
    # llnull = const
    const_in = np.ones_like(score)
    X = np.stack([const_in, score], axis=1)

    try:
        logit = sm.Logit(phe, X).fit()

        print('regression')
        # print(logit.summary())
        print('coef:\n', logit.params)

        llf = logit.llf
        llnull = logit.llnull
    except np.linalg.LinAlgError as e:
        print(e)
        llf, llnull = np.nan, np.nan
    except BaseException as e:
        print('another error: ', e)
        # include 'PerfectSeparationError
        llf, llnull = np.nan, np.nan

    return llf, llnull



def calc_likelihood_metric(phe,score,score_cov):

    llmodel, _ = calc_likelihood_regression_score(score,phe)
    llnull, _ = calc_likelihood_regression_score(score_cov, phe)

    return llmodel, llnull 


def nagelkerke_r2(phe, score, score_cov):
    # no randomness
    np.random.seed(51)

    llmodel, llnull = calc_likelihood_metric(phe,score,score_cov)

    nsample = len(phe)
    # nagel = (1 - (lnull / lmodel)**(2 / nsample)) / (1 - lnull**(2 / nsample))
    nagel = (1 - np.exp((llnull - llmodel) * (2 / nsample))) / (1 - np.exp(llnull * (2 / nsample)))
    print('nagel', nagel)
    
    return nagel """