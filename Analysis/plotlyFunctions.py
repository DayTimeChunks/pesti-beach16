import plotly.graph_objs as go

from scipy import stats
import numpy as np
import statsmodels.api as sm


def get_xy_line(df, obs, delta=False):
    if delta:
        min_x = min(df['SIM'])
        min_y = min(df[obs])
        min_plot = min(min_x, min_y)
        max_plot = max(max(df['SIM']), max(df[obs]))
        x_test = [min_plot - 1, max_plot + 1]
        y_test = x_test

        trace_bck = go.Scatter(
            x=x_test,
            y=y_test,
            name='y = x',
            mode='lines',
            marker=dict(
                color='grey',
                opacity=0.9
            )
        )
    else:
        max_x = max(df['SIM'])
        max_y = max(df[obs])
        max_plot = max(max_x, max_y)
        min_plot = 0

        x_test = [min_plot, max_plot]
        y_test = x_test

        trace_bck = go.Scatter(
            x=x_test,
            y=y_test,
            name='y = x',
            mode='lines',
            marker=dict(
                color='grey',
                opacity=0.9
            )
        )
    return trace_bck


def get_scatter(df, obs, sd, name=None):
    labels = df['IDcal'].tolist()
    sim = df['SIM'].tolist()
    obs = df[obs].tolist()
    err = df[sd].tolist()

    trace = go.Scatter(
        x=sim,
        y=obs,
        name=name,
        error_y=dict(
            type='data',
            symmetric=True,
            array=err,
            # arrayminus=errMin
        ),
        mode='markers+text',
        text=labels,
        textposition='right',  # 'bottom center',
        textfont=dict(
            # family='sans serif',
            size=7,
            # color='#ff7f0e'
        )
    )

    return trace


def get_layout(df, var, folders, title='Title', delta=False, x_jitter=0, y_jitter=0):
    # obs = df[var].tolist()
    # sim = df['SIM'].tolist()

    r = ''
    n = ''
    for folder in range(1, folders + 1):
        df1 = df.loc[df['IDsim'] == ('Sim-' + str(folder))]

        # R-squared
        obs = df1[var].tolist()
        sim = df1['SIM'].tolist()
        r_squared = str('R' + str(folder) + ' = ' + str(round(get_r_squared(obs, sim), 2)) + '\n \n')
        r += r_squared

        # Nash
        mean = df1[var].mean()
        # Diff sim vs. obs
        dfn = df1.assign(diff_sim=lambda row: (row['SIM'] - row[var]) ** 2)
        err_sim = dfn['diff_sim'].sum()
        # Variance
        dfn = dfn.assign(diff_obs=lambda row: (row[var] - mean) ** 2)
        err_obs = dfn['diff_obs'].sum()
        error = err_sim / err_obs
        nash = 1 - error
        if not delta:  # Log only for concentrations
            lnmean = np.log(df1[var]).mean()
            # Log Diff sim vs. obs
            dfn = dfn.assign(lndiff_sim=lambda row: (np.log(row['SIM']) - np.log(row[var])) ** 2)
            err_lnsim = dfn['lndiff_sim'].sum()
            # Log variance
            dfn = dfn.assign(lndiff_obs=lambda row: (np.log(row[var]) - lnmean) ** 2)
            err_lnobs = dfn['lndiff_obs'].sum()
            error += err_lnsim / err_lnobs
            nash = 1 - 0.5*error

        n += str('N' + str(folder) + ' = ' + str(round(nash, 2)) + '\n \n')

    if delta:
        x = max(df['SIM']) + 0.7
        y = min(df[var]) - 0.7
        y2 = y + 0.7
    else:
        x = max(df[var]) - 0.1 * max(df[var]) + x_jitter
        y = min(df['SIM']) + 0.1 * max(df['SIM'])
        y2 = y + 1 + y_jitter

    layout = go.Layout(
        title=title,
        xaxis=dict(
            title='Simulated',
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        ),
        yaxis=dict(
            title='Observed',
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        ),
        showlegend=True,
        annotations=[
            dict(
                x=x,
                y=y,
                # xref='x',
                # yref='y',
                text=r,
                showarrow=False
                # arrowhead=7,
                # xshift = 10,
                # ax=0,
                # ay=-40
                # opacity = 0.7
            ),
            dict(
                x=x,
                y=y2,
                # xref='x',
                # yref='y',
                text=n,
                showarrow=False
                # arrowhead=7,
                # xshift = 10,
                # ax=0,
                # ay=-40
                # opacity = 0.7
            )
        ]
    )

    return layout


def get_r_squared(obs, sim):
    # Add y-intercept to regression data
    obs.insert(0, 0)
    sim.insert(0, 0)

    x = sm.add_constant(np.array(obs), prepend=True)  # Declares first element as y-intercept
    model = sm.OLS(np.array(sim), x)
    results = model.fit()

    return results.rsquared
