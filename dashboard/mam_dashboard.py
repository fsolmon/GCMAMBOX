import pathlib
import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.colors as mpc
import numpy as np
import xarray as xr
import subprocess
import os
import tempfile
import uuid
import shutil

# ============================================================================
# Log-normal distribution functions (from original script)
# ============================================================================

def lognormal_distribution(Dg, sigma_g, N, x_min=0.001, x_max=100, num_points=100):
    """Generate x, y values for a log-normal distribution"""
    mu = np.log(Dg)
    sigma = np.log(sigma_g)
    x = np.logspace(np.log10(x_min), np.log10(x_max), num_points)
    y = (N / (np.sqrt(2 * np.pi) * sigma)) * \
        np.exp(-((np.log(x) - mu)**2) / (2 * sigma**2))
    return x, y

def lognormal_volume_distribution(Dg, sigma_g, V_total, x_min=0.001, x_max=100., num_points=100):
    """Generate x, y values for a log-normal volume distribution"""
    mu_number = np.log(Dg)
    sigma = np.log(sigma_g)
    Dg_volume = Dg * np.exp(3 * sigma**2)
    mu_volume = np.log(Dg_volume)
    x = np.logspace(np.log10(x_min), np.log10(x_max), num_points)
    y = (V_total / (np.sqrt(2 * np.pi) * sigma)) * \
        np.exp(-((np.log(x) - mu_volume)**2) / (2 * sigma**2))
    return x, y

# ============================================================================
# Constants and color schemes
# ============================================================================

gascolor = {
    'so2_gas': 'green',
    'h2so4_gas': 'lightgreen',
    'hno3_gas': 'blue',
    'nh3_gas': 'orange',
    'soag_gas': 'indigo',
    'hcl_gas': 'cyan',
}

def coltr(colo,t):
    if (colo == 'toto'):
      out = 'rgba'+ str(mpc.to_rgba(colo)[:3]+(t,))
    else:
      out = colo
    return out
tr = 0.8
specolor = {
 'bc_aer':  coltr('black', tr),
 'pom_aer': coltr('magenta', tr),
 'dst_aer': coltr('yellow', tr),
 'co3_aer': coltr('burlywood', tr),
'ca_aer':  coltr('pink', tr),
 'ncl_aer': coltr('grey', tr),
  'cl_aer':  coltr('cyan', tr),
 'so4_aer': coltr('green',tr),
    'no3_aer': coltr('blue',tr),
    'nh4_aer': coltr('orange', tr),
    'soa_aer': coltr('indigo',tr),
'wat_aer': coltr('lightblue', 0.6),

    }

modname = {0: 'Accumulation', 1: 'Aitken', 2: 'Coarse', 3: 'Primary Carbon'}

sigmag = np.array([1.800, 1.600, 1.800, 1.600])

# ============================================================================
# Utility functions - ONLY MODIFIED write_namelist and run_fortran_model
# ============================================================================

def write_namelist(params, namelist_path='namelist'):
    """Write namelist file from parameters dictionary"""
    namelist_content = f"""&time_input
mam_dt         = {params['mam_dt']},
mam_nstep      = {params['mam_nstep']},
/
&cntl_input
mdo_mambox     = {params['mdo_mambox']},
mdo_coldstart  = {params['mdo_coldstart']},
mdo_gaschem    = {params['mdo_gaschem']},
mdo_gasaerexch = {params['mdo_gasaerexch']},
mdo_rename     = {params['mdo_rename']},
mdo_newnuc     = {params['mdo_newnuc']},
mdo_coag       = {params['mdo_coag']},
/
&met_input
mtmin           = {params['mtmin']},
mtmax           = {params['mtmax']},
mrhmin          = {params['mrhmin']},
mrhmax         = {params['mrhmax']},
press          = {params['press']},
/
&chem_input
numc = {params['numc'][0]},{params['numc'][1]},{params['numc'][2]},{params['numc'][3]}
mfso4 = {params['mfso4'][0]}, {params['mfso4'][1]}, {params['mfso4'][2]}, {params['mfso4'][3]}
mfpom = {params['mfpom'][0]}, {params['mfpom'][1]}, {params['mfpom'][2]}, {params['mfpom'][3]}
mfsoa = {params['mfsoa'][0]}, {params['mfsoa'][1]}, {params['mfsoa'][2]}, {params['mfsoa'][3]}
mfbc  = {params['mfbc'][0]}, {params['mfbc'][1]}, {params['mfbc'][2]}, {params['mfbc'][3]}
mfdst = {params['mfdst'][0]}, {params['mfdst'][1]}, {params['mfdst'][2]}, {params['mfdst'][3]}
mfncl = {params['mfncl'][0]}, {params['mfncl'][1]}, {params['mfncl'][2]}, {params['mfncl'][3]}
mfno3 = {params['mfno3'][0]}, {params['mfno3'][1]}, {params['mfno3'][2]}, {params['mfno3'][3]}
mfnh4 = {params['mfnh4'][0]}, {params['mfnh4'][1]}, {params['mfnh4'][2]}, {params['mfnh4'][3]}
mfco3 = {params['mfco3'][0]}, {params['mfco3'][1]}, {params['mfco3'][2]}, {params['mfco3'][3]}
mfca  = {params['mfca'][0]}, {params['mfca'][1]}, {params['mfca'][2]}, {params['mfca'][3]}
mfcl  = {params['mfcl'][0]}, {params['mfcl'][1]}, {params['mfcl'][2]}, {params['mfcl'][3]}
qso2          = {params['qso2']},
qh2so4        = {params['qh2so4']},
qsoag         = {params['qsoag']},
qhno3         = {params['qhno3']},
qnh3          = {params['qnh3']},
qhcl          = {params['qhcl']}
/
"""
    with open(namelist_path, 'w') as f:
        f.write(namelist_content)

def run_fortran_model(work_dir='.'):
    """Run the Fortran gcmambox executable"""

    BASE_DIR = pathlib.Path(__file__).resolve().parent.parent
    executable = BASE_DIR / "modbin" / "gcmambox"
    try:
        result = subprocess.run([str(executable)],
                              capture_output=True,
                              text=True,
                              timeout=30,
                              cwd=work_dir))
        if result.returncode != 0:
            return False, f"Error: {result.stderr}"
        return True, "Simulation completed successfully ----- Click play for results !"
    except subprocess.TimeoutExpired:
        return False, "Simulation timed out (>30s)"
    except FileNotFoundError:
        return False, "gcmambox executable not found"
    except Exception as e:
        return False, f"Error running simulation: {str(e)}"

def create_animated_figure(nc_file='mam_output.nc', dt = 1200):
    """Create the animated plotly figure from NetCDF output"""
    try:
        ds = xr.open_dataset(nc_file)
        print('open', nc_file)
    except FileNotFoundError:
        return go.Figure().add_annotation(
            text="No simulation data available. Run simulation first.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )

    frames = []
    for s in ds.nsteps:
        df = ds.sel(nsteps=s)
        data_list = []
        ntrac1 = 0
        # Gas concentrations
        for g in list(gascolor.keys()):
            conc = ds[g].sel(nsteps=slice(0, int(s+1)))
            trace = go.Scatter(
                x=conc.nsteps.values,
                y=conc.values,
                name=g,
                line_color=gascolor[g]
            )
            data_list.append(trace)
            ntrac1 = ntrac1 +1
        # Volume distribution (stacked)
        fil = 'tozeroy'
        ycum = 0
        ntrac2 = ntrac1
        for c in list(specolor.keys()):
            ymodcum = 0.
            for m in df.mode.values:
                Dg = df['Dgn_mode'].sel(mode=m).values * 1.E6
                sigma_g = sigmag[m]
                mass = df[c].sel(mode=m).values
                x_vol, y_vol = lognormal_volume_distribution(Dg, sigma_g, mass)
                ymodcum = ymodcum + y_vol
            ycum = ycum + ymodcum
            trace = go.Scatter(
                x=x_vol, y=ycum, fill=fil,
                mode='none', fillcolor=specolor[c],
                opacity=0.3, name=c
            )
            data_list.append(trace)
            ntrac2 = ntrac2+1
            fil = 'tonexty'

        # Total mass
        trace = go.Scatter(
            x=x_vol, y=ycum, mode='lines',
            name='Total Mass', line_color='grey'
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1
# Text
        trace = go.Scatter(
        x=[0.005] ,
        y=[np.max(ycum)*0.6] ,
#        y = [5.E-9],
        mode='text',  # Show both markers and text
        text=['RH = %s'%round(float(df['relhum']),2)],
        name = "   ",
        textposition='top left',
#        marker=dict(size=10, color='blue'),
        textfont=dict(size=16, color='blue')
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1

        trace = go.Scatter(
        x=[0.005] ,
        y=[np.max(ycum)*0.8] ,
#        y = [12.E-9],
        mode='text',  # Show both markers and text
        text=['TEMP = %s'%round(float(df['temp']),2)],
        textposition='top left',
        name= "  ",
        #        marker=dict(size=10, color='blue'),
        textfont=dict(size=16, color='red')
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1

        trace = go.Scatter(
        x=[50] ,
        y=[np.max(ycum)*0.8] ,
#        y = [12.E-9],
        mode='text',  # Show both markers and text
        text=[],
        name = "  ",
        textposition='top left',
#        marker=dict(size=10, color='blue'),
        textfont=dict(size=16, color='red')
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1

# Number distribution
        ymodcum = 0
        ntrac3 = ntrac2
        for m in df.mode.values:
            Dg = df['Dgn_mode'].sel(mode=m).values * 1.E6
            sigma_g = sigmag[m]
            num = df['num_aer'].sel(mode=m).values
            x_num, y_num = lognormal_distribution(Dg, sigma_g, num)
            trace = go.Scatter(
                x=x_num, y=y_num, fill='tozeroy',
                mode='none', name=modname[int(m)]
            )
            data_list.append(trace)
            ymodcum = ymodcum + y_num
            ntrac3 = ntrac3+1
        trace = go.Scatter(
            x=x_num, y=ymodcum, mode='lines',
            name='Total number', line_color='black'
        )
        data_list.append(trace)
        ntrac3 = ntrac3+1


        frames += [go.Frame(data=data_list, name='step%s' % s)]

    # Create figure with subplots
    fig = go.Figure().set_subplots(3, 1,vertical_spacing=0.07)
                                   #insets=[{'cell':(2,1), 'l':0.5, 'b':0.5 , 'w' : 0.1, 'h' : 0.1}])
#    fig = make_subplots(rows=3, cols=1,
#                    specs=[[{'type': 'xy'}], [{'type': 'xy'}], [{'type': 'xy'}]],
#                    vertical_spacing=0.07,
#                    insets=[dict(cell=(1,1), l=0.55, b= 0.43, w = 0.1, h = 0.1),
#                            dict(cell=(2,1), l=0.5, h=0.65, b=0.1,  type='polar')])
    fig.add_traces(frames[0].data[0:ntrac1], 1, 1)
    fig.add_traces(frames[0].data[ntrac1:ntrac2], 2, 1)
    fig.add_traces(frames[0].data[ntrac2:ntrac3], 3, 1)



    # Update axes
    fig.update_xaxes(range=[-1, len(ds.nsteps)], title='time steps', row=1, col=1)
    fig.update_yaxes(type='log', range=[-1., 2], title='gas conc. (ppb)', row=1, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=2, col=1)
    fig.update_yaxes(range=[0, np.max(ycum)], title='aer. dm/dlogD (kg.kg-1)', row=2, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=3, col=1)
    fig.update_yaxes(type='log', range=[5, 11], title='aer. dN/dlogD (#.kg-1)', row=3, col=1)

    fig.frames = frames

# Animation controls
    fr_duration = 500
    sliders = [{
        "pad": {"b": 10, "t": 10},
        "len": 0.9,
        "x": 0.1,
        "y": 1.0,  # Move to top
        "yanchor": "bottom",
        "currentvalue":{"prefix": "Time (mn) : "},
        "steps": [
            {
                "args": [[f.name], {
                    "frame": {"duration": fr_duration},
                    "mode": "immediate",
                    "fromcurrent": True,
                    "transition": {"duration": fr_duration, "easing": "linear"},
                }],
#                "label": f"tstep{k+1}",add time label in mn
                "label": (k)*dt / 60.,
                "method": "animate",
            }
            for k, f in enumerate(fig.frames)
        ],
    }]

    fig.update_layout(
        sliders=sliders,
        legend={'traceorder': 'normal','font': {'size': 18} },
        height=1200,
        margin=dict(t=80, b=40, l=60, r=40),  # Add top margin for controls
        updatemenus=[{
            "buttons": [
                {
                    "args": [None, {
                        "frame": {"duration": fr_duration},
                        "mode": "immediate",
                        "fromcurrent": True,
                        "transition": {"duration": fr_duration, "easing": "linear"},
                    }],
                    "label": "&#9654;",
                    "method": "animate",
                },
                {
                    "args": [[None], {
                        "frame": {"duration": fr_duration},
                        "mode": "immediate",
                        "fromcurrent": True,
                        "transition": {"duration": fr_duration, "easing": "linear"},
                    }],
                    "label": "&#9724;",
                    "method": "animate",
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 10},
            "type": "buttons",
            "x": 0.0,
            "xanchor": "left",
            "y": 1.02,  # Move to top
            "yanchor": "bottom",
        }]
    )

    return fig

# ============================================================================
# Helper function for creating mass fraction inputs
# ============================================================================

def create_mass_fraction_input(species, mode, default_value):
    """Create a compact input field for mass fractions"""
    return html.Div([
        html.Label(f"{species}:", style={'fontSize': '11px', 'width': '35px', 'display': 'inline-block'}),
        dcc.Input(
            id=f'mf{species.lower()}-{mode}',
            type='number',
            value=default_value,
            min=0,
            max=1,
            step=0.01,
            style={'width': 'calc(100% - 40px)', 'fontSize': '11px'}
        ),
    ], style={'marginBottom': '2px'})

# ============================================================================
# Dash App Layout
# ============================================================================

app = dash.Dash(__name__)

app.layout = html.Div([

#    html.Div([
#           html.Img(src='assets/logo.png', style={'height': '130px', 'width': '600px'}),],
#           style={'textAlign': 'center', 'marginBottom': '10px', 'padding': '10px'}),

html.Div([
    html.Div(style={'flex': '1'}),  # Spacer
    html.Img(src='assets/logo.png', style={'height': '130px', 'width': '419px'}),
    html.Div([
        html.A('📄 Doc !',
               href='assets/doc.pdf',
               target='_blank',
               style={'padding': '16px 16px',
                      'backgroundColor': '#0066cc', 'color': 'white', 'textDecoration': 'none',
                      'borderRadius': '20px', 'fontSize': '16px'})
    ], style={'flex': '1', 'textAlign': 'right', 'display': 'flex', 'alignItems': 'center', 'justifyContent': 'flex-end'})
], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '10px', 'padding': '10px'}),


html.Div([
        # Left panel - Controls (scrollable)
        html.Div([
            html.Div([
                html.H3("Simulation Controls", style={'marginTop': '0'}),

                # Run button at the top
                html.Button('Run Simulation', id='run-button', n_clicks=0,
                           style={'width': '100%', 'height': '40px', 'fontSize': '16px',
                                 'backgroundColor': '#4CAF50', 'color': 'white',
                                 'border': 'none', 'cursor': 'pointer', 'borderRadius': '5px', 'marginBottom': '10px'}),

                html.Div(id='status-message', style={'marginTop': '5px', 'marginBottom': '10px', 'fontWeight': 'bold', 'fontSize': '12px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Time parameters
                html.H4("Time Parameters", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("Time step (s):", style={'fontSize': '12px'}),
                dcc.Input(id='mam_dt', type='number', value=1200, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Number of steps:", style={'fontSize': '12px'}),
                dcc.Input(id='mam_nstep', type='number', value=100, style={'width': '100%', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Process controls
                html.H4("Process Controls", style={'fontSize': '16px', 'marginTop': '10px'}),
                dcc.Checklist(
                    id='process_controls',
                    options=[
                        {'label': ' Gas Chemistry', 'value': 'gaschem'},
                        {'label': ' Gas-Aerosol Exchange', 'value': 'gasaerexch'},
                        {'label': ' Rename', 'value': 'rename'},
                        {'label': ' New Nucleation', 'value': 'newnuc'},
                        {'label': ' Coagulation', 'value': 'coag'},
                    ],
                    value=['gaschem', 'gasaerexch', 'rename', 'newnuc', 'coag'],
                    style={'fontSize': '12px'}
                ),

                html.Hr(style={'margin': '10px 0'}),

                # Meteorological parameters
                html.H4("Meteorological Parameters", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("Temperature Initial (K):", style={'fontSize': '12px'}),
                dcc.Input(id='mtmin', type='number', value=275, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Temperature Final (K):", style={'fontSize': '12px'}),
                dcc.Input(id='mtmax', type='number', value=275, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Pressure (Pa):", style={'fontSize': '12px'}),
                dcc.Input(id='press', type='number', value=1.e5, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Relative Humidity Initial:", style={'fontSize': '12px'}),
                dcc.Input(id='mrhmin', type='number', value=0.5, step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Relative Humidity Final:", style={'fontSize': '12px'}),
                dcc.Input(id='mrhmax', type='number', value=0.5, step=0.01, style={'width': '100%', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Gas concentrations
                html.H4("Gas Concentrations (ppb)", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("SO2:", style={'fontSize': '12px'}),
                dcc.Input(id='qso2', type='number', value=30., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("H2SO4:", style={'fontSize': '12px'}),
                dcc.Input(id='qh2so4', type='number', value=0., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("HNO3:", style={'fontSize': '12px'}),
                dcc.Input(id='qhno3', type='number', value=30., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("NH3:", style={'fontSize': '12px'}),
                dcc.Input(id='qnh3', type='number', value=5., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("SOAG:", style={'fontSize': '12px'}),
                dcc.Input(id='qsoag', type='number', value=0., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("HCl:", style={'fontSize': '12px'}),
                dcc.Input(id='qhcl', type='number', value=0., style={'width': '100%', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Number concentrations by mode
                html.H4("Number Concentrations (#/cm³)", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Div([
                    html.Div([
                        html.Label("Accum:", style={'fontSize': '11px', 'width': '50px', 'display': 'inline-block'}),
                        dcc.Input(id='numc-0', type='number', value=0.e5, style={'width': 'calc(100% - 55px)', 'fontSize': '11px'}),
                    ], style={'marginBottom': '3px'}),
                    html.Div([
                        html.Label("Aitken:", style={'fontSize': '11px', 'width': '50px', 'display': 'inline-block'}),
                        dcc.Input(id='numc-1', type='number', value=0.e5, style={'width': 'calc(100% - 55px)', 'fontSize': '11px'}),
                    ], style={'marginBottom': '3px'}),
                    html.Div([
                        html.Label("Coarse:", style={'fontSize': '11px', 'width': '50px', 'display': 'inline-block'}),
                        dcc.Input(id='numc-2', type='number', value=0.e7, style={'width': 'calc(100% - 55px)', 'fontSize': '11px'}),
                    ], style={'marginBottom': '3px'}),
                    html.Div([
                        html.Label("Prim C:", style={'fontSize': '11px', 'width': '50px', 'display': 'inline-block'}),
                        dcc.Input(id='numc-3', type='number', value=0.e5, style={'width': 'calc(100% - 55px)', 'fontSize': '11px'}),
                    ], style={'marginBottom': '10px'}),
                ]),

                html.Hr(style={'margin': '10px 0'}),

                # Mass fractions
                html.H4("Mass Fractions by Mode", style={'fontSize': '16px', 'marginTop': '10px'}),

                # Mode 0 - Accumulation
                html.Details([
                    html.Summary("Mode 0 - Accumulation", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 0, 0.2),
                        create_mass_fraction_input("POM", 0, 0.1),
                        create_mass_fraction_input("SOA", 0, 0.1),
                        create_mass_fraction_input("BC", 0, 0.1),
                        create_mass_fraction_input("DST", 0, 0.1),
                        create_mass_fraction_input("NCL", 0, 0.15),
                        create_mass_fraction_input("NO3", 0, 0.0),
                        create_mass_fraction_input("NH4", 0, 0.0),
                        create_mass_fraction_input("CO3", 0, 0.0),
                        create_mass_fraction_input("CA", 0, 0.0),
                        create_mass_fraction_input("CL", 0, 0.15),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 1 - Aitken
                html.Details([
                    html.Summary("Mode 1 - Aitken", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 1, 0.8),
                        create_mass_fraction_input("POM", 1, 0.0),
                        create_mass_fraction_input("SOA", 1, 0.2),
                        create_mass_fraction_input("BC", 1, 0.0),
                        create_mass_fraction_input("DST", 1, 0.0),
                        create_mass_fraction_input("NCL", 1, 0.0),
                        create_mass_fraction_input("NO3", 1, 0.0),
                        create_mass_fraction_input("NH4", 1, 0.0),
                        create_mass_fraction_input("CO3", 1, 0.0),
                        create_mass_fraction_input("CA", 1, 0.0),
                        create_mass_fraction_input("CL", 1, 0.0),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 2 - Coarse
                html.Details([
                    html.Summary("Mode 2 - Coarse", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 2, 0.0),
                        create_mass_fraction_input("POM", 2, 0.0),
                        create_mass_fraction_input("SOA", 2, 0.0),
                        create_mass_fraction_input("BC", 2, 0.0),
                        create_mass_fraction_input("DST", 2, 0.0),
                        create_mass_fraction_input("NCL", 2, 0.4),
                        create_mass_fraction_input("NO3", 2, 0.0),
                        create_mass_fraction_input("NH4", 2, 0.0),
                        create_mass_fraction_input("CO3", 2, 0.0),
                        create_mass_fraction_input("CA", 2, 0.0),
                        create_mass_fraction_input("CL", 2, 0.6),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 3 - Primary Carbon
                html.Details([
                    html.Summary("Mode 3 - Primary Carbon", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 3, 0.0),
                        create_mass_fraction_input("POM", 3, 0.5),
                        create_mass_fraction_input("SOA", 3, 0.0),
                        create_mass_fraction_input("BC", 3, 0.5),
                        create_mass_fraction_input("DST", 3, 0.0),
                        create_mass_fraction_input("NCL", 3, 0.0),
                        create_mass_fraction_input("NO3", 3, 0.0),
                        create_mass_fraction_input("NH4", 3, 0.0),
                        create_mass_fraction_input("CO3", 3, 0.0),
                        create_mass_fraction_input("CA", 3, 0.0),
                        create_mass_fraction_input("CL", 3, 0.0),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '10px'}),

            ], style={'padding': '15px', 'overflowY': 'scroll', 'height': 'calc(100vh - 130px)'})
        ], style={'width': '15%', 'backgroundColor': '#f5f5f5', 'boxSizing': 'border-box'}),

        # Right panel - Visualization
        html.Div([
            dcc.Loading(
                id="loading",
                type="default",
                children=dcc.Graph(id='animation-plot', style={'height': 'calc(100vh - 160px)'})
            ),
        ], style={'width': '85%', 'padding': '10px', 'boxSizing': 'border-box', 'backgroundColor': '#e6f2ff'}),

    ], style={'display': 'flex', 'height': 'calc(100vh - 160px)'}),
], style={'fontFamily': 'Arial, sans-serif', 'margin': '0', 'padding': '0', 'backgroundColor': '#e6f2ff'})

# ============================================================================
# Callbacks - ONLY MODIFIED THE SIMULATION EXECUTION PART
# ============================================================================

@app.callback(
    [Output('animation-plot', 'figure'),
     Output('status-message', 'children')],
    [Input('run-button', 'n_clicks')],
    [State('mam_dt', 'value'),
     State('mam_nstep', 'value'),
     State('process_controls', 'value'),
     State('mtmin', 'value'),
     State('mtmax', 'value'),
     State('press', 'value'),
     State('mrhmin', 'value'),
     State('mrhmax', 'value'),
     State('qso2', 'value'),
     State('qh2so4', 'value'),
     State('qhno3', 'value'),
     State('qnh3', 'value'),
     State('qsoag', 'value'),
     State('qhcl', 'value'),
     # Number concentrations
     State('numc-0', 'value'),
     State('numc-1', 'value'),
     State('numc-2', 'value'),
     State('numc-3', 'value'),
     # Mass fractions - Mode 0
     State('mfso4-0', 'value'),
     State('mfpom-0', 'value'),
     State('mfsoa-0', 'value'),
     State('mfbc-0', 'value'),
     State('mfdst-0', 'value'),
     State('mfncl-0', 'value'),
     State('mfno3-0', 'value'),
     State('mfnh4-0', 'value'),
     State('mfco3-0', 'value'),
     State('mfca-0', 'value'),
     State('mfcl-0', 'value'),
     # Mass fractions - Mode 1
     State('mfso4-1', 'value'),
     State('mfpom-1', 'value'),
     State('mfsoa-1', 'value'),
     State('mfbc-1', 'value'),
     State('mfdst-1', 'value'),
     State('mfncl-1', 'value'),
     State('mfno3-1', 'value'),
     State('mfnh4-1', 'value'),
     State('mfco3-1', 'value'),
     State('mfca-1', 'value'),
     State('mfcl-1', 'value'),
     # Mass fractions - Mode 2
     State('mfso4-2', 'value'),
     State('mfpom-2', 'value'),
     State('mfsoa-2', 'value'),
     State('mfbc-2', 'value'),
     State('mfdst-2', 'value'),
     State('mfncl-2', 'value'),
     State('mfno3-2', 'value'),
     State('mfnh4-2', 'value'),
     State('mfco3-2', 'value'),
     State('mfca-2', 'value'),
     State('mfcl-2', 'value'),
     # Mass fractions - Mode 3
     State('mfso4-3', 'value'),
     State('mfpom-3', 'value'),
     State('mfsoa-3', 'value'),
     State('mfbc-3', 'value'),
     State('mfdst-3', 'value'),
     State('mfncl-3', 'value'),
     State('mfno3-3', 'value'),
     State('mfnh4-3', 'value'),
     State('mfco3-3', 'value'),
     State('mfca-3', 'value'),
     State('mfcl-3', 'value'),
    ],
    prevent_initial_call=False
)
def run_simulation(n_clicks, mam_dt, mam_nstep, processes, mtmin,mtmax, press, mrhmin,mrhmax,
                   qso2, qh2so4, qhno3, qnh3, qsoag, qhcl,
                   numc0, numc1, numc2, numc3,
                   # Mode 0 mass fractions
                   mfso4_0, mfpom_0, mfsoa_0, mfbc_0, mfdst_0, mfncl_0,
                   mfno3_0, mfnh4_0, mfco3_0, mfca_0, mfcl_0,
                   # Mode 1 mass fractions
                   mfso4_1, mfpom_1, mfsoa_1, mfbc_1, mfdst_1, mfncl_1,
                   mfno3_1, mfnh4_1, mfco3_1, mfca_1, mfcl_1,
                   # Mode 2 mass fractions
                   mfso4_2, mfpom_2, mfsoa_2, mfbc_2, mfdst_2, mfncl_2,
                   mfno3_2, mfnh4_2, mfco3_2, mfca_2, mfcl_2,
                   # Mode 3 mass fractions
                   mfso4_3, mfpom_3, mfsoa_3, mfbc_3, mfdst_3, mfncl_3,
                   mfno3_3, mfnh4_3, mfco3_3, mfca_3, mfcl_3):
    """Run simulation when button is clicked and update the plot"""

    # For initial load, just create figure from existing data
    if n_clicks == 0:
        if os.path.exists('mam_output.nc'):
            fig = create_animated_figure()
            return fig, "Showing previous simulation results"
        else:
            empty_fig = go.Figure().add_annotation(
                text="Click 'Run Simulation' to start",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=24)
            )
            return empty_fig, "Ready to run simulation"

    # Create unique temporary directory for this simulation
    session_id = str(uuid.uuid4())
    temp_dir = tempfile.mkdtemp(prefix=f'mam_sim_{session_id}_')

    try:
        # Build parameters dictionary
        params = {
            'mam_dt': mam_dt,
            'mam_nstep': mam_nstep,
            'mdo_gaschem': 1 if 'gaschem' in processes else 0,
            'mdo_gasaerexch': 1 if 'gasaerexch' in processes else 0,
            'mdo_rename': 1 if 'rename' in processes else 0,
            'mdo_newnuc': 1 if 'newnuc' in processes else 0,
            'mdo_coag': 1 if 'coag' in processes else 0,
            'mdo_mambox': 1,
            'mdo_coldstart': 1,
            'mtmin': mtmin,
            'mtmax': mtmax,
            'mrhmin':mrhmin,
            'mrhmax':mrhmax,
            'press': press,
            'numc': [numc0, numc1, numc2, numc3],
            'mfso4': [mfso4_0, mfso4_1, mfso4_2, mfso4_3],
            'mfpom': [mfpom_0, mfpom_1, mfpom_2, mfpom_3],
            'mfsoa': [mfsoa_0, mfsoa_1, mfsoa_2, mfsoa_3],
            'mfbc': [mfbc_0, mfbc_1, mfbc_2, mfbc_3],
            'mfdst': [mfdst_0, mfdst_1, mfdst_2, mfdst_3],
            'mfncl': [mfncl_0, mfncl_1, mfncl_2, mfncl_3],
            'mfno3': [mfno3_0, mfno3_1, mfno3_2, mfno3_3],
            'mfnh4': [mfnh4_0, mfnh4_1, mfnh4_2, mfnh4_3],
            'mfco3': [mfco3_0, mfco3_1, mfco3_2, mfco3_3],
            'mfca': [mfca_0, mfca_1, mfca_2, mfca_3],
            'mfcl': [mfcl_0, mfcl_1, mfcl_2, mfcl_3],
            'qso2': qso2,
            'qh2so4': qh2so4,
            'qsoag': qsoag,
            'qhno3': qhno3,
            'qnh3': qnh3,
            'qhcl': qhcl,
        }

        # Copy executable to temp directory
        if os.path.exists('./gcmambox'):
            shutil.copy('./gcmambox', temp_dir)
        else:
            return (go.Figure().add_annotation(
                text="gcmambox executable not found",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=20, color='red')
            ), "❌ Executable not found")

        # Write namelist file in temp directory
        namelist_path = os.path.join(temp_dir, 'namelist')
        write_namelist(params, namelist_path)

        # Run Fortran model in temp directory
        success, message = run_fortran_model(temp_dir)

        if not success:
            error_fig = go.Figure().add_annotation(
                text=f"Simulation failed: {message}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=20, color='red')
            )
            return error_fig, f"❌ {message}"

        # Copy output file back to main directory with unique session name
        output_nc = os.path.join(temp_dir, 'mam_output.nc')
        if os.path.exists(output_nc):
            session_output = f'mam_output_{session_id}.nc'
            shutil.copy(output_nc, session_output)
            print('mam_dt',mam_dt)
            # Create figure from session-specific results
            fig = create_animated_figure(session_output, dt=mam_dt)

            # Clean up session-specific output file after creating figure
            try:
                os.remove(session_output)
            except Exception as e:
                print(f"Warning: Could not remove {session_output}: {e}")

            return fig, "✅ " + message
        else:
            return (go.Figure().add_annotation(
                text="No output file generated",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=20, color='red')
            ), "❌ No output generated")

    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            print(f"Warning: Could not clean up temp directory {temp_dir}: {e}")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)
