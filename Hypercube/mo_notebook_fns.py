import pandas as pd
import numpy as np
from SALib.analyze import morris as moa

# Return an array of Si's (dictionary: mu, sigma, param_names, etc)
def getSiList(ouput_keys, problem, vector_test, grid_jump, p, runs,
              nash=True, console=False):
    morris_results = []
    for indx in range(len(ouput_keys)):
        var_name = ouput_keys[indx]
        # One Y-array per variable, len(Y) = no. of runs
        Y = getTestOutput(var_name, runs, nash=nash)
        Si = moa.analyze(problem, vector_test, Y,
                         num_resamples=1000, print_to_console=console,
                         grid_jump=grid_jump, num_levels=p)
        morris_results.append(Si)
    return morris_results

# Return a Y-array (one variable only, with length == runs)
def getTestOutput(var, runs, nash=True):
    Y = np.zeros([runs])
    for run in range(runs):
        path = str(run + 1) + "/" + var + ".tss"
        res = np.loadtxt(path, float, skiprows=4, usecols=[1])
        if nash:
            Y[run] = res[-1] # Take the last Nash value
        else:
            Y[run] = res.sum() # Take the cumulative
    return Y


# Save Morris as a JSON file
def saveMorris(ouput_keys, SiList, filename):
    import json
    mListSi = dict()
    for i in range(len(ouput_keys)):
        var_name = ouput_keys[i]
        Si = SiList[i]
        mSi = dict()
        mSi["mu"] = Si["mu"].tolist()
        mSi["mu_star"] = Si["mu_star"].tolist()
        mSi['mu_star_conf'] = Si['mu_star_conf']
        mSi['sigma'] = Si['sigma'].tolist()
        mSi['names'] = Si['names']
        mListSi[var_name] = mSi

    with open(filename, 'w') as f:
        json.dump(mListSi, f)

# Used in get_morris_df()
def getMorris(keys, morris_results):
    """
    :param keys:
    :param morris_results:
    :return: A dictionary of morris evaluations
    """
    mListSi = dict()
    for i in range(len(keys)):
        var_name = keys[i]
        Si = morris_results[i]
        mSi = dict()
        mSi["mu"] = Si["mu"].tolist()
        mSi["mu_star"] = Si["mu_star"].tolist()
        mSi['mu_star_conf'] = Si['mu_star_conf']
        mSi['sigma'] = Si['sigma'].tolist()
        mSi['names'] = Si['names']
        mListSi[var_name] = mSi

    return mListSi


def get_morris_df(tss_names, results, norm=True):
    """
    :param tss_names: TSS list to convert to dataframes
    :param results:
    :param norm:
    :return: A concatenated dataframe of all tss passed in
    """
    mDict = getMorris(tss_names, results)
    df_array = []
    for name in range(len(tss_names)):
        df = pd.DataFrame({"parameter": mDict[tss_names[name]]["names"],
                      "mu_star": mDict[tss_names[name]]["mu_star"],
                       "sigma": mDict[tss_names[name]]["sigma"]
                      })
        df['name'] = tss_names[name]
        # Normalize
        if norm:
            df['mu_star'] = df['mu_star']/max(df['mu_star'])
            df['sigma'] = df['sigma']/max(df['sigma'])
        df_array.append(df)
    mDf = pd.concat(df_array)
    return mDf