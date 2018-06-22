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
            mode='lines',
            marker=dict(
                color='grey',
                opacity=0.9
            )
        )
    return trace_bck


def get_scatter(df, obs, sd):

    labels = df['IDcal'].tolist()
    sim = df['SIM'].tolist()
    obs = df[obs].tolist()
    err = df[sd].tolist()

    trace = go.Scatter(
        x=sim,
        y=obs,
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


def get_layout(df, var, title='Title', x_jitter=0, y_jitter=0):

    obs = df[var].tolist()
    sim = df['SIM'].tolist()

    x = max(df[var]) - (0.1 + y_jitter) * max(df[var])
    y = min(df['SIM']) + (0.1 + x_jitter) * max(df['SIM'])

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
        showlegend=False,
        annotations=[
            dict(
                x=x,
                y=y,
                # xref='x',
                # yref='y',
                text=str('R-sq. = ' + str(round(get_r_squared(obs, sim), 2))),
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
