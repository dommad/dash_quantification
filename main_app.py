"""Dash app for label-free quantification of proteomics data"""
from dash import Dash, dcc, html, dash_table, Input, Output
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import quantify as qu

app = Dash(__name__)

def get_res(files, imput=""):
    """combine results from multiple files"""
    option = 1
    if imput == "Average":
        option = 0
    results = qu.Quantify().combine_replicates(files)
    qu.Quantify().imput_missing(results, option=option)
    return results

# load data from the files
files1 = ("interact-f10.ipro.prot.xml", "interact-f10.ipro.prot.xml")
files2 = ("interact-f07.ipro.prot.xml", "interact-f07.ipro.prot.xml")
f10 = get_res(files1)
f07 = get_res(files2)
f10_ave = get_res(files1, imput="Average")
f07_ave = get_res(files2, imput="Average")

app = Dash(__name__)

app.layout = html.Div(
    style={
        'display': 'flex',
        'margin-left': '20px'
        },
        children = [
    html.Div(style={'display': 'inline-block', 'width': '50%'},
    children=[
    html.H1("Quantification of Proteomics Data"),
    html.Br(),
    html.H2("Quantification method"),
    dcc.RadioItems(
        ['NSAF', 'SI'],
        'SI',
        id='select-quant-method'
        ),
    html.Br(),
    html.H2("Select FDR threshold"),
    dcc.Slider(
        0, 0.05, step=None,
        id='select-fdr',
        marks=dict(
            zip(
                np.linspace(0, 0.05, 11),
                [str(i) for i in np.linspace(0, 0.05, 11)]
                )
                ),
        value=0.01
    ),
    html.Br(),
    html.H2("Select log2(fold change) (absolute)"),
    dcc.Slider(id='select-fold-change',
            min=1,
            max=5,
            marks={i: str(i) for i in np.linspace(1, 5, 9)},
            value=1,
        ),
    html.Br(),
    html.H2("Select imputation method for missing values"),
    dcc.Dropdown(
        ['Average', 'Minimum in replicate'],
        'Minimum in replicate',
        id='select-imputation'
        )]),

    html.Div(
        style={
            'display': 'inline-block',
            'width': '50%',
            'margin-left': '20px',
            'margin-right': '20px',
            'margin-top': '20px'
            },
        children=[
            html.Div(
                id='volcano-plot'
                    ),
            html.Div(
                style={
                    'margin-left': '20px',
                    'margin-right': '20px',
                    'margin-top': '20px'},
                id='data-table'
                    )
])]
)

@app.callback(
[Output('volcano-plot', 'children'),
Output('data-table', 'children')],
Input('select-fdr', 'value'),
Input('select-fold-change', 'value'),
Input('select-quant-method', 'value'),
Input('select-imputation', 'value'))
def update_volcano(fdr, fold, q_mode, imput):
    """update the volcano plot"""
    if imput == "Average":
        new_d = qu.Quantify().welch_test(f10_ave, f07_ave, mode=q_mode)
    else:
        new_d = qu.Quantify().welch_test(f10, f07, mode=q_mode)

    sel_data = qu.Quantify().get_stats(new_d, fdr, fold)

    new_fig = go.Figure(data=[
                        go.Scattergl(
                        y=[-np.log(el[1][0]) for el in sel_data],
                        x=[el[1][1] for el in sel_data],
                        mode='markers',
                        text=[x[0] for x in sel_data],
                        )
                        ],
                        layout={
                            "xaxis": {'title': "log2(fold change)"},
                            "yaxis": {'title': "-log(p-value)"}
                            }
                        )

    dash_t = get_table(sel_data)
    return dcc.Graph(figure=new_fig), dash_t


def get_table(sel_data):
    """generate dash table from the quantification data provided"""
    tab_df = pd.DataFrame.from_dict(
        dict(sel_data),
        orient='index',
        columns=['p-value', 'log2(fold change)']
        )
    tab_df['protein'] = tab_df.index
    tab_df.reset_index(inplace=True, drop=True)

    dash_t = dash_table.DataTable(
        id='table-proteins',
        columns=[
                {"name": i, "id": i, "deletable": True, "selectable": True} for i in tab_df.columns
                ],
        data=tab_df.to_dict('records'),
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        row_selectable="multi",
        row_deletable=True,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 10,
        style_cell={
                    'backgroundColor':'rgba(223, 234, 185, 0.23)',
                    'opacity':0.95,
                    'color': 'black'
                    }
        )
    return dash_t



if __name__ == '__main__':
    app.run_server(debug=True)
