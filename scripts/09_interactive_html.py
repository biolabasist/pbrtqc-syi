"""Interactive full-year control charts (self-contained HTML).

Usage: python scripts/09_interactive_html.py <output directory>
"""

import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join('..', 'output')
html_dir = os.path.join(out, 'interactive')

for k in ('creatinine', 'sodium'):
    d = pd.read_csv(os.path.join(html_dir, f'chart_{k}.csv'))
    label, config, ucl, lcl = open(os.path.join(html_dir, f'chart_{k}.txt'), encoding='utf-8').read().split('\n')[:4]
    t = pd.to_datetime(d['time']).to_numpy()
    y = d['stat'].to_numpy(np.float32)
    a = d['alarm'].to_numpy(bool)
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=t, y=y, mode='lines', name='Control statistic',
                               line=dict(color='#08519c', width=1),
                               hovertemplate='%{x|%Y-%m-%d %H:%M}<br>Statistic: %{y:.2f}%<extra></extra>'))
    fig.add_trace(go.Scattergl(x=t[a], y=y[a], mode='markers', name='Out of control',
                               marker=dict(color='#E8112D', size=3), hoverinfo='skip'))
    for v in (float(ucl), float(lcl)):
        fig.add_hline(y=v, line=dict(color='#E8112D', dash='dash', width=1))
    fig.update_layout(title=f'{label}, {config}, 2023: control statistic as % deviation from the in-control mean',
                      xaxis=dict(title='Result time', rangeslider=dict(visible=True)),
                      yaxis=dict(title='% from in-control mean'), template='plotly_white')
    fig.write_html(os.path.join(html_dir, f'Supplementary_interactive_chart_{k}.html'),
                   include_plotlyjs=True, full_html=True)
