"""Dash app for label-free quantification of proteomics data"""
import argparse
from dash import Dash, dcc, html, dash_table, Input, Output, State
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import quantify as qu



def execute_main(c_files, t_files):

    controls_min = qu.generate_master_dict(c_files, impute_method="minimum")
    treatments_min = qu.generate_master_dict(t_files, impute_method="minimum")
    controls_ave = qu.generate_master_dict(c_files, impute_method="average")
    treatments_ave = qu.generate_master_dict(t_files, impute_method="average")

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
            marks={int(i) if i % 1 == 0 else i: '{}'.format(i) for i in np.linspace(0, 0.05, 11)},
            value=0.01
        ),
        html.Br(),
        html.H2("Select log2(fold change) (absolute)"),
                dcc.Slider(
                           0.0, 5.0, step=None, id='select-fold-change',
                marks={int(i) if i % 1 == 0 else i: '{}'.format(i) for i in np.linspace(0, 5, 11)},
                value=1.5
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
                        ),

    ])]
    )


    @app.callback(
    [Output('volcano-plot', 'children'),
    Output('data-table', 'children')],
    Input('select-fdr', 'value'),
    Input('select-fold-change', 'value'),
    Input('select-quant-method', 'value'),
    Input('select-imputation', 'value'))
    def update_volcano(fdr_lim, fold_lim, q_mode, imput):
        """update the volcano plot"""
        qmodes = {'SI': 'sin', 'NSAF': 'nsaf'}

        if imput == "Average":
            new_d = qu.welch_test(controls_ave, treatments_ave, qmodes[q_mode])
        else:
            new_d = qu.welch_test(controls_min, treatments_min, qmodes[q_mode])

        sel_data = qu.get_stats(new_d, fdr_lim, fold_lim)

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
                        },
            export_format='xlsx',
            export_headers='display',
            )
        return dash_t
    
    app.run_server(debug=True)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Quantification of shotgun proteomics data')
    parser.add_argument('-c',  '--control', nargs='*', required=True,
                        type=str, help='prot.xml file(s) output by StPeter containing data for control replicates')
    parser.add_argument('-t', '--treatment', nargs='*', required=True,
                        type=str, help="prot.xml file(s) output by StPeter containing data for treatment replicates")
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    c_files = args.control
    t_files = args.treatment
    output_name = args.output
    execute_main(c_files, t_files)