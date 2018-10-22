import numpy as np
import pandas as pd
import os

def concat_periods(measure, periods, path):
    """
    Concat B2L (coupled) periods (tss files)
    :param measure:
    :param periods:
    :param path:
    :return: None, saves to CSV
    """
    if measure == 'Q':
        # Simulated
        filename = "resW_accQ_m3.tss"
        col = 'Q Coup.'
    elif measure == 'LF':
        filename = "resW_outLatflow_m3.tss"
        col = 'LF Coup.'
    elif measure == 'ROFF':
        filename = "resW_accRunoff_m3.tss"
        col = 'ROFF Coup.'
    elif measure == 'ADR':
        filename = "resW_accDrain_m3.tss"
        col = 'ADR Coup.'
    elif measure == 'BF':
        filename = "resW_accBaseflow_m3.tss"
        col = 'BF Coup.'
    elif measure == 'Bal_W':
        filename = "resW_global_waterMB.tss"
        col = 'BalW Coup.'

    dsets = []
    for p in range(periods):
        sim_path = path + "\\output_coup\\" + str(p) + "\\"
        sim = pd.read_table(sim_path + filename,
                            skiprows=4, delim_whitespace=True,
                            names=['Jdays', col],
                            header=None)
        if p + 1 < periods:
            sim = sim.iloc[:-2, :]
        dsets.append(sim)

    sim = pd.concat(dsets)
    # Save file
    all_path = path + "\\output_coup\\all\\"
    if not os.path.exists(all_path):
        os.mkdir(all_path)
    sim.to_csv(all_path + filename, sep='\t')


def get_obs(name):
    """
    :param name:
    :return: obs data as pandas data frame
    """
    obs_path = './observations/' + name + ".tss"
    return pd.read_table(obs_path)


def get_dataframe(sim_path, measure):
    """
    :param sim_path:
    :param measure:
    :return: data frame for get_kge()
    """
    # Observed
    name_obs = measure.lower() + '_cal'
    obs = get_obs(name_obs)
    col = 'sim'
    if measure == 'Q_out':
        # Simulated
        filename = "resW_accQ_m3.tss"
        col = 'Q Sim.'
    elif measure == 'CONC_out':
        # Simulated
        filename = "resM_oCONC_ugL.tss"
    elif measure == 'd13C_out':
        # Simulated
        filename = "resM_outISO_d13C.tss"
    elif measure == 'LDS_out':
        # Simulated
        filename = "resM_EXP_light_g.tss"
    else:
        print("No appropriate measure selected for the outlet!")
        return

    sim = pd.read_table(sim_path + filename,
                        skiprows=4, delim_whitespace=True,
                        names=['Jdays', col],
                        header=None)
    # Merge
    match = pd.merge(obs, sim, how='inner', on='Jdays')
    return match


def get_kge(df, col_obs, col_sim):
    """
    :param df:
    :param col: name of observation column (e.g. Qm3)
    :return:
    """
    # KGE, @Gupta2009
    # Linear correlation, alpha and beta
    r = np.corrcoef(df[col_obs], df[col_sim])[1, 0]
    alpha = np.std(df[col_sim]) / np.std(df[col_obs])
    beta = np.mean(df[col_sim]) / np.mean(df[col_obs])
    kge = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    return kge


