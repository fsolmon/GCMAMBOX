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

# Resolve the executable and optics files relative to this script, not the launch directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GCMAMBOX_EXE = os.path.join(SCRIPT_DIR, 'gcmambox')
ASSETS_DIR = os.path.join(SCRIPT_DIR, 'assets')
DOC_FILE = os.path.join(SCRIPT_DIR, '..', 'doc', 'doc_gcmambox.md')
PRESET_DIR = os.path.join(ASSETS_DIR, 'presets')

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
 'na_aer': coltr('grey', tr),
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
accom_coef_h2so4 = {params['accom_coef_h2so4']},
accom_coef_hno3  = {params['accom_coef_hno3']},
accom_coef_hcl   = {params['accom_coef_hcl']},
accom_coef_nh3   = {params['accom_coef_nh3']},
/
&met_input
mtmin           = {params['mtmin']},
mtmax           = {params['mtmax']},
mrhmin          = {params['mrhmin']},
mrhmax         = {params['mrhmax']},
press          = {params['press']},
/
&chem_input
masstot = {params['masstot'][0]},{params['masstot'][1]},{params['masstot'][2]},{params['masstot'][3]}
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
&emis_input
h2so4_prod    = {params['h2so4_prod']},
hno3_prod     = {params['hno3_prod']},
nh3_prod      = {params['nh3_prod']},
hcl_prod      = {params['hcl_prod']},
soag_prod     = {params['soag_prod']},
/
"""
    with open(namelist_path, 'w') as f:
        f.write(namelist_content)

def run_fortran_model(work_dir='.'):
    """Run the Fortran gcmambox executable"""
    try:
        result = subprocess.run(['./gcmambox'],
                              capture_output=True,
                              text=True,
                              timeout=30,
                              cwd=work_dir)
        if result.returncode != 0:
            return False, f"Error: {result.stderr}"
        return True, "Simulation completed successfully ----- Click play for results !"
    except subprocess.TimeoutExpired:
        return False, "Simulation timed out (>30s)"
    except FileNotFoundError:
        return False, "gcmambox executable not found"
    except Exception as e:
        return False, f"Error running simulation: {str(e)}"

def get_optical_input_files(base_dir='./'):
    """Get list of optical property NetCDF files needed by MAM"""
    # List of required optical input files
    optical_files = [
        # Species optical properties
        'sulfate_rrtmg_c080918.nc',
        'bcpho_rrtmg_c100508.nc',
        'ocphi_rrtmg_c100508.nc',
        'ocpho_rrtmg_c130709.nc',
        'dust_aeronet_rrtmg_c141106.nc',
        'ssam_rrtmg_c100508.nc',
        # Modal optical properties
        'mam4_mode1_rrtmg_aeronetdust_c141106.nc',
        'mam4_mode2_rrtmg_aitkendust_c141106.nc',
        'mam4_mode3_rrtmg_aeronetdust_c141106.nc',
        'mam4_mode4_rrtmg_c130628.nc',
    ]
    
    # Check which files exist; return full paths of existing ones and names of missing ones
    existing_files = []
    missing_files = []
    for fname in optical_files:
        fpath = os.path.join(base_dir, fname)
        if os.path.exists(fpath):
            existing_files.append(fpath)
        else:
            print(f"Error: Optical file not found: {fpath}")
            missing_files.append(fname)

    return existing_files, missing_files

# ============================================================================
# Preset cases: namelist files in assets/presets/ (also runnable directly by gcmambox)
# Header lines "! title: ..." and "! description: ..." are shown in the dashboard.
# ============================================================================

MF_SPECIES = ['so4', 'pom', 'soa', 'bc', 'dst', 'ncl', 'no3', 'nh4', 'co3', 'ca', 'cl']
PROCESS_FLAGS = [('mdo_gasaerexch', 'gasaerexch'), ('mdo_rename', 'rename'),
                 ('mdo_newnuc', 'newnuc'), ('mdo_coag', 'coag')]
PRESET_SCALARS = ['mam_dt', 'mam_nstep', 'mtmin', 'mtmax', 'press', 'mrhmin', 'mrhmax',
                  'qh2so4', 'qhno3', 'qnh3', 'qsoag', 'qhcl',
                  'accom_coef_h2so4', 'accom_coef_hno3', 'accom_coef_hcl', 'accom_coef_nh3',
                  'hno3_prod', 'nh3_prod', 'soag_prod', 'hcl_prod']
# dashboard input ids filled by a preset
PRESET_FIELDS = (['process_controls', 'h2so4_prod'] + PRESET_SCALARS
                 + [f'masstot-{m}' for m in range(4)]
                 + [f'mf{sp}-{m}' for sp in MF_SPECIES for m in range(4)])

def read_preset(path):
    """Parse a preset namelist file -> (title, description, {name: float or list})"""
    title, description, values = os.path.basename(path), '', {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith('! title:'):
                title = line.split(':', 1)[1].strip()
            elif line.lower().startswith('! description:'):
                description = line.split(':', 1)[1].strip()
            line = line.split('!')[0].strip()
            if '=' not in line or line.startswith('&'):
                continue
            key, val = line.split('=', 1)
            nums = [float(v) for v in val.replace(',', ' ').split()]
            values[key.strip().lower()] = nums[0] if len(nums) == 1 else nums
    return title, description, values

def list_presets():
    """Dropdown options for the preset files, sorted by file name"""
    if not os.path.isdir(PRESET_DIR):
        return []
    return [{'label': read_preset(os.path.join(PRESET_DIR, f))[0], 'value': f}
            for f in sorted(os.listdir(PRESET_DIR)) if f.endswith('.nml')]

def preset_to_inputs(values):
    """Map preset namelist values to dashboard input values (ordered as PRESET_FIELDS)"""
    out = {}
    if any(flag in values for flag, _ in PROCESS_FLAGS):
        out['process_controls'] = [name for flag, name in PROCESS_FLAGS if values.get(flag, 0) > 0]
    if 'h2so4_prod' in values:
        out['h2so4_prod'] = values['h2so4_prod'] / 1000.   # ppt/h in the namelist, ppb/h in the UI
    for key in PRESET_SCALARS:
        if key in values:
            out[key] = int(values[key]) if key in ('mam_dt', 'mam_nstep') else values[key]
    for m, v in enumerate(values.get('masstot', [])):
        out[f'masstot-{m}'] = v
    for sp in MF_SPECIES:
        for m, v in enumerate(values.get(f'mf{sp}', [])):
            out[f'mf{sp}-{m}'] = v
    return [out.get(field, dash.no_update) for field in PRESET_FIELDS]

def create_animated_figure(nc_file='mam_output.nc', dt=1200, press=1.e5):
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
    maxyaod = ds['aod_mode'].max(dim='nsteps').max(dim='mode').values * 1.4

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
                name=f'{g} ({float(df[g].values):.3g})',
                line_color=gascolor[g],
                legend='legend'
            )
            data_list.append(trace)
            ntrac1 = ntrac1 +1

        # RH and TEMP gauges — top-right of subplot 1 (gas concentrations)
        rh_value   = float(df['relhum'] * 100)
        temp_value = float(df['temp'])

        x_rh_pos   = len(ds.nsteps) * 0.78
        x_temp_pos = len(ds.nsteps) * 0.92
        y_rh_pos   = 10**1.15  # ~14 ppb on log scale, clear of the panel top
        y_temp_pos = 10**1.15

        # RH colour: white → blue
        if rh_value > 75:
            rh_outer = '#0047AB'; rh_inner = '#1E90FF'; text_color_rh = 'white'
        elif rh_value > 50:
            rh_outer = '#4169E1'; rh_inner = '#6495ED'; text_color_rh = 'white'
        elif rh_value > 25:
            rh_outer = '#87CEEB'; rh_inner = '#B0E0E6'; text_color_rh = 'black'
        else:
            rh_outer = '#E6F3FF'; rh_inner = '#FFFFFF';  text_color_rh = 'black'

        # TEMP colour: white → red
        if temp_value > 315:
            temp_outer = '#8B0000'; temp_inner = '#DC143C'; text_color_temp = 'black'
        elif temp_value > 295:
            temp_outer = '#FF4500'; temp_inner = '#FF6347'; text_color_temp = 'black'
        elif temp_value > 273:
            temp_outer = '#FFA07A'; temp_inner = '#FFB6C1'; text_color_temp = 'black'
        else:
            temp_outer = '#FFE6E6'; temp_inner = '#FFFFFF';  text_color_temp = 'black'

        trace_rh_bg = go.Scatter(
            x=[x_rh_pos, x_rh_pos, x_rh_pos],
            y=[y_rh_pos, y_rh_pos, y_rh_pos],
            mode='markers',
            marker=dict(
                size=[80, 75, 65],
                color=[rh_outer, 'white', rh_inner],
                opacity=[0.9, 1.0, 0.95],
                line=dict(color=rh_outer, width=[3, 2, 2])
            ),
            showlegend=False, hoverinfo='skip', name=''
        )
        data_list.append(trace_rh_bg)
        ntrac1 = ntrac1 + 1

        trace_rh_text = go.Scatter(
            x=[x_rh_pos], y=[y_rh_pos],
            mode='text',
            text=[f'<b>RH</b><br>{rh_value:.0f}<br><b>%</b>'],
            textfont=dict(size=12, color=text_color_rh, family='Arial Black'),
            textposition='middle center',
            showlegend=False, hoverinfo='skip', name=''
        )
        data_list.append(trace_rh_text)
        ntrac1 = ntrac1 + 1

        trace_temp_bg = go.Scatter(
            x=[x_temp_pos, x_temp_pos, x_temp_pos],
            y=[y_temp_pos, y_temp_pos, y_temp_pos],
            mode='markers',
            marker=dict(
                size=[80, 75, 65],
                color=[temp_outer, 'white', temp_inner],
                opacity=[0.9, 1.0, 0.95],
                line=dict(color=temp_outer, width=[3, 2, 2])
            ),
            showlegend=False, hoverinfo='skip', name=''
        )
        data_list.append(trace_temp_bg)
        ntrac1 = ntrac1 + 1

        trace_temp_text = go.Scatter(
            x=[x_temp_pos], y=[y_temp_pos],
            mode='text',
            text=[f'<b>TEMP</b><br>{temp_value:.0f}<br><b>K</b>'],
            textfont=dict(size=12, color=text_color_temp, family='Arial Black'),
            textposition='middle center',
            showlegend=False, hoverinfo='skip', name=''
        )
        data_list.append(trace_temp_text)
        ntrac1 = ntrac1 + 1

        # data is already in µg/m³ (conversion done in driver.F90)
        conv = 1.0

        # Volume distribution (stacked)
        fil = 'tozeroy'
        ycum = 0
        ntrac2 = ntrac1

        # Pre-compute integrated masses in µg/m³: total and per-mode
        species_mode_mass = {}
        species_total_mass = {}
        for c in list(specolor.keys()):
            m_masses = [float(df[c].sel(mode=m).values) * conv for m in df.mode.values]
            species_mode_mass[c] = m_masses
            species_total_mass[c] = sum(m_masses)
        grand_total = sum(species_total_mass.values())
        # total mass per mode (sum over all species), keyed by mode coordinate value
        mode_total_mass = {int(m): sum(species_mode_mass[c][i] for c in specolor.keys())
                           for i, m in enumerate(df.mode.values)}

        for c in list(specolor.keys()):
            ymodcum = 0.
            for m in df.mode.values:
                Dg = df['Dgn_mode'].sel(mode=m).values * 1.E6
                sigma_g = sigmag[m]
                mass = float(df[c].sel(mode=m).values) * conv
                x_vol, y_vol = lognormal_volume_distribution(Dg, sigma_g, mass)
                ymodcum = ymodcum + y_vol
            ycum = ycum + ymodcum

            sp_total = species_total_mass[c]
            m_masses = species_mode_mass[c]
            significant = grand_total > 0 and sp_total / grand_total > 0.01
            trace_name = f'{c} ({sp_total:.2f})' if significant else c

            # customdata columns: [total, mode0, mode1, ...] in µg/m³
            cdata = np.tile([sp_total] + m_masses, (len(x_vol), 1))
            mode_hover = ''.join(
                f'{modname[int(m)]}: %{{customdata[{i+1}]:.2f}} µg/m³<br>'
                for i, m in enumerate(df.mode.values)
            )
            htemplate = (
                f'<b>{c}</b><br>'
                f'Total: %{{customdata[0]:.2f}} µg/m³<br>'
                + mode_hover
                + '<extra></extra>'
            )

            trace = go.Scatter(
                x=x_vol, y=ycum, fill=fil,
                mode='none', fillcolor=specolor[c],
                opacity=0.3, name=trace_name,
                legend='legend2',
                customdata=cdata, hovertemplate=htemplate
            )
            data_list.append(trace)
            ntrac2 = ntrac2+1
            fil = 'tonexty'

        # Total mass envelope — per-mode totals in hover
        mode_totals = [
            sum(species_mode_mass[c][i] for c in specolor.keys())
            for i in range(len(df.mode.values))
        ]
        total_cdata = np.tile([grand_total] + mode_totals, (len(x_vol), 1))
        total_htemplate = (
            f'<b>Total Mass</b>: %{{customdata[0]:.2f}} µg/m³<br>'
            + ''.join(
                f'{modname[int(m)]}: %{{customdata[{i+1}]:.2f}} µg/m³<br>'
                for i, m in enumerate(df.mode.values)
            )
            + '<extra></extra>'
        )
        trace = go.Scatter(
            x=x_vol, y=ycum, mode='lines',
            name=f'Total Mass ({grand_total:.2f})', line_color='grey',
            legend='legend2',
            customdata=total_cdata, hovertemplate=total_htemplate
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1

# Number distribution, converted from #/kg-air to #/cm3 with the air density at this step
        kg_to_cm3 = press / (287.05 * float(df['temp'])) * 1.e-6
        ymodcum = 0
        ntrac3 = ntrac2
        for m in df.mode.values:
            Dg = df['Dgn_mode'].sel(mode=m).values * 1.E6
            sigma_g = sigmag[m]
            num = df['num_aer'].sel(mode=m).values * kg_to_cm3
            x_num, y_num = lognormal_distribution(Dg, sigma_g, num)
            trace = go.Scatter(
                x=x_num, y=y_num, fill='tozeroy',
                mode='none', name=f'{modname[int(m)]} ({float(num):.3g})',
                legend='legend3'
            )
            data_list.append(trace)
            ymodcum = ymodcum + y_num
            ntrac3 = ntrac3+1
        trace = go.Scatter(
            x=x_num, y=ymodcum, mode='lines',
            name=f"Total number ({float(df['num_aer'].sum().values) * kg_to_cm3:.3g})", line_color='black',
            legend='legend3'
        )
        data_list.append(trace)
        ntrac3 = ntrac3+1
              
#  Modal AOD Bar Chart
        ntrac4 = ntrac3
        
        # Check if AOD data exists in dataset
        if 'aod_mode' in ds.variables:
            # Get wet diameters for bar X positions (convert to micrometers)
            x_positions = ['Aitken', 'PC','Accum','Coarse' ]
            #for m in df.mode.values:
            #    Dg_wet = df['Dgn_mode'].sel(mode=m).values * 1.E6
            #    x_positions.append(Dg_wet)
            #    # Calculate bar widths for constant visual width on log scale
            #log_width = 0.25  # Width in log10 units (adjust to taste)
            #bar_widths = []
            #for x_center in x_positions:
            #    log_center = np.log10(x_center)
            #    x_left = 10**(log_center - log_width/2)
            #    x_right = 10**(log_center + log_width/2)
            #    bar_widths.append(x_right - x_left)
            
            # Define species and colors for stacking
            aod_species = ['aod_sulfate', 'aod_bc', 'aod_pom','aod_mom', 'aod_soa', 'aod_dust', 'aod_seasalt','aod_no3',  'aod_nh4']
            species_colors = {
                'aod_sulfate': 'green',
                'aod_bc': 'black',
                'aod_pom': 'magenta',
                'aod_soa': 'indigo',
                'aod_dust': 'yellow',
                'aod_seasalt': 'grey',
                'aod_mom': 'brown',      # Marine organic matter
                'aod_no3': 'blue',       
                'aod_nh4': 'orange',     
            }
            
            # Create stacked bars for each species
            for species in aod_species:
                if species in ds.variables:
                    y_values = []
             #       for m in df.mode.values:
                    for m in [1,3,0,2]:
                       aod_val = df[species].sel(mode=m).values
                       y_values.append(float(aod_val))
                    
                    # Create bar trace
                    trace = go.Bar(
                        x=x_positions,
                        y=y_values,
                        name=species, #.replace('aod_', '').upper(),
                        marker_color=species_colors[species],
                        opacity=0.75,
                        legend='legend4'
                        #width=bar_widths,  # Array of widths for constant log appear
                    )
                    data_list.append(trace)
                    ntrac4 = ntrac4 + 1
            
            # Add total SSA arkers at top of stacks
            total_aod_values = []
            total_ssa_values = []
            total_g_values   = []
            total_hyg_values = []
            total_lwabs_values = []
            total_pH_values = []
            mode_mass_aod = []
            ytextpos_values = []
            ytextpos_ssa_values = []
            ytextpos_g_values = []
            ytextpos_hyg_values = []
            ytextpos_lwabs_values = []
            ytextpos_pH_values = []
            #for m in df.mode.values:
            print(df['hygro_aer'].values)
            for m in [1,3,0,2]:
                #if 'aod_mode' in ds.variables:
                    total_aod = df['aod_mode'].sel(mode=m).values 
                    total_aod_values.append(float(total_aod))
                    total_ssa =  df['ssa_mode'].sel(mode=m).values
                    total_ssa_values.append(float(total_ssa))
                    total_g =  df['asm_mode'].sel(mode=m).values
                    total_g_values.append(float(total_g))
                    total_hyg =  df['hygro_aer'].sel(mode=m).values
                    total_hyg_values.append(float(total_hyg))
                    if 'aod_lw_mode' in ds.variables:
                        total_lwabs = df['aod_lw_mode'].sel(mode=m).values
                        total_lwabs_values.append(float(total_lwabs))
                    ytextpos = 0.85 * maxyaod  
                    ytextpos_values.append(float(ytextpos)) 
                    ytextpos = 0.80 * maxyaod
                    ytextpos_ssa_values.append(float(ytextpos))
                    ytextpos = 0.75 * maxyaod
                    ytextpos_g_values.append(float(ytextpos))
                    ytextpos = 0.91 * maxyaod
                    ytextpos_hyg_values.append(float(ytextpos))
                    ytextpos = 0.70 * maxyaod
                    ytextpos_lwabs_values.append(float(ytextpos))
                    if 'pH_aer' in ds.variables:
                        # pH_aer is already actual pH (-log10 of the linearly-averaged
                        # H+ concentration), computed in driver.F90; NaN means no aerosol
                        # water / no H+ diagnosed for that mode.
                        total_pH_values.append(float(df['pH_aer'].sel(mode=m).values))
                    mode_mass_aod.append(mode_total_mass.get(int(m), 0.0))
                    ytextpos = 0.65 * maxyaod
                    ytextpos_pH_values.append(float(ytextpos))
                #else:
                #    total_aod_values.append(0.0)
                #    mode_ssa_values.append(0.0)
            # ssa as scatter plot on top
            trace = go.Scatter(
                x=x_positions,
               # y=total_aod_values ,switch to fixed position 
                y= ytextpos_hyg_values , 
                mode='text',
                name='mode hygrosc.',
                marker=dict(size=14, color='black', symbol='diamond'),
             #   text=['ssa = '+ f'{v:.3f}' for v in mode_ssa_values],
                text=['hygros. = '+ f'{v:.3f}' for v in total_hyg_values],
                textposition='top center',
                textfont=dict(size=16, color='blue'),
                showlegend=False
            )
            data_list.append(trace)
            ntrac4 = ntrac4 + 1

            trace = go.Scatter(
                x=x_positions,
               # y=total_aod_values ,switch to fixed position
                y= ytextpos_values ,
                mode='text',
                name='Mode AOD',
                marker=dict(size=14, color='black', symbol='diamond'),
             #   text=['ssa = '+ f'{v:.3f}' for v in mode_ssa_values],
                text=['aod = '+ f'{v:.3f}' for v in total_aod_values],
                textposition='top center',
                textfont=dict(size=16, color='black'),
                showlegend=False
            )
            data_list.append(trace)
            ntrac4 = ntrac4 + 1

            trace = go.Scatter(
                x=x_positions,
               # y=total_aod_values ,switch to fixed position 
                y= ytextpos_ssa_values ,
                mode='text',
                name='Mode ssa',
                marker=dict(size=14, color='black', symbol='diamond'),
             #   text=['ssa = '+ f'{v:.3f}' for v in mode_ssa_values],
                text=['ssa = '+ f'{v:.3f}' for v in total_ssa_values],
                textposition='top center',
                textfont=dict(size=16, color='black'),
                showlegend=False
            )
            data_list.append(trace)
            ntrac4 = ntrac4 + 1

            trace = go.Scatter(
                x=x_positions,
               # y=total_aod_values ,switch to fixed position 
                y= ytextpos_g_values ,
                mode='text',
                name='Mode asy',
                marker=dict(size=14, color='black', symbol='diamond'),
             #   text=['ssa = '+ f'{v:.3f}' for v in mode_ssa_values],
                text=['asy = '+ f'{v:.3f}' for v in total_g_values],
                textposition='top center',
                textfont=dict(size=16, color='black'),
                showlegend=False
            )
            data_list.append(trace)
            ntrac4 = ntrac4 + 1

            if 'aod_lw_mode' in ds.variables and total_lwabs_values:
                trace = go.Scatter(
                    x=x_positions,
                    y=ytextpos_lwabs_values,
                    mode='text',
                    name='abs LW',
                    text=['LW aaod = ' + f'{v:.3f}' for v in total_lwabs_values],
                    textposition='top center',
                    textfont=dict(size=16, color='red'),
                    showlegend=False
                )
                data_list.append(trace)
                ntrac4 = ntrac4 + 1

            if 'pH_aer' in ds.variables and total_pH_values:
                trace = go.Scatter(
                    x=x_positions,
                    y=ytextpos_pH_values,
                    mode='text',
                    name='pH',
                    text=['pH = ' + (f'{v:.2f}' if not np.isnan(v) and mode_mass_aod[i] > 0 else 'N/A')
                          for i, v in enumerate(total_pH_values)],
                    textposition='top center',
                    textfont=dict(size=16, color='purple'),
                    showlegend=False
                )
                data_list.append(trace)
                ntrac4 = ntrac4 + 1

        else:
            # Placeholder if no AOD data
            trace = go.Scatter(
                x=[0.1], y=[0.1],
                mode='text',
                text=['AOD data not available'],
                textfont=dict(size=14, color='red'),
                showlegend=False
            )
            data_list.append(trace)
            ntrac4 = ntrac4 + 1


        frames += [go.Frame(data=data_list, name='step%s' % s)]

    # Create figure with subplots
    fig = go.Figure().set_subplots(4, 1, vertical_spacing=0.07, row_heights=[0.22, 0.26, 0.22, 0.30],
                                   subplot_titles=('<b>Gas concentrations</b>',
                                                   '<b>Aerosol mass size distribution</b>',
                                                   '<b>Aerosol number size distribution</b>',
                                                   '<b>Mode properties</b>'))

    fig.add_traces(frames[0].data[0:ntrac1], 1, 1)
    fig.add_traces(frames[0].data[ntrac1:ntrac2], 2, 1)
    fig.add_traces(frames[0].data[ntrac2:ntrac3], 3, 1)
    fig.add_traces(frames[0].data[ntrac3:ntrac4], 4, 1)



    # Update axes
    fig.update_xaxes(range=[-1, len(ds.nsteps)], title='time steps', row=1, col=1)
    fig.update_yaxes(type='log', range=[-1., 2], title='gas conc. (ppb)', row=1, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=2, col=1)
    fig.update_yaxes(range=[0, np.max(ycum)], title='aer. dm/dlogD (µg.m⁻³)', row=2, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=3, col=1)
    fig.update_yaxes(type='log', range=[-1, 6], title='aer. dN/dlogD (#.cm⁻³)', row=3, col=1)

    
#    fig.update_xaxes(type="log", range=[-3, 2], title='Wet Diameter (µm)', row=4, col=1)
    fig.update_xaxes(title='MAM mode', tickfont=dict(size=16), row=4, col=1)
    fig.update_yaxes(range=[0, maxyaod ], title='wet AOD @ 550nm / DZ =1 km', row=4, col=1)

    fig.frames = frames

# Animation controls
    fr_duration = 500
    sliders = [{
        "pad": {"b": 10, "t": 10},
        "len": 0.9,
        "x": 0.1,
        "y": 1.025,  # above the first subplot title
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
        barmode='stack',  # Enable stacking for AOD bars
        height=1600,
        margin=dict(t=120, b=40, l=60, r=40),  # top margin for the play/slider controls
        # Legends aligned with the top of their subplot (row domains from set_subplots)
        # Legend 1 - Gas concentrations (top subplot)
        legend=dict(
            x=1.02, y=fig.layout.yaxis.domain[1], xanchor='left', yanchor='top',
            font=dict(size=14),
            traceorder='normal',
            title=dict(text='<b>Gas Concentration (ppb)</b>', font=dict(size=14))
        ),
        # Legend 2 - Mass distribution
        legend2=dict(
            x=1.02, y=fig.layout.yaxis2.domain[1], xanchor='left', yanchor='top',
            font=dict(size=14),
            traceorder='normal',
            title=dict(text='<b>Mass Concentration (µg/m³)</b>', font=dict(size=14))
        ),
        # Legend 3 - Number distribution
        legend3=dict(
            x=1.02, y=fig.layout.yaxis3.domain[1], xanchor='left', yanchor='top',
            font=dict(size=14),
            traceorder='normal',
            title=dict(text='<b>Number Concentration (#/cm³)</b>', font=dict(size=14))
        ),
        # Legend 4 - AOD (bottom subplot)
        legend4=dict(
            x=1.02, y=fig.layout.yaxis4.domain[1], xanchor='left', yanchor='top',
            font=dict(size=14),
            traceorder='normal'
        ),
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
            "y": 1.045,  # above the first subplot title
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

@app.server.route('/doc')
def serve_doc():
    """Render doc/doc_gcmambox.md as an HTML page (linked from the Doc button)"""
    import html as htmllib
    try:
        with open(DOC_FILE, encoding='utf-8') as f:
            text = f.read()
    except OSError:
        return f"Documentation file not found: {DOC_FILE}", 404
    try:
        import markdown
        body = markdown.markdown(text, extensions=['tables', 'fenced_code'])
    except ImportError:
        body = f"<pre>{htmllib.escape(text)}</pre>"
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>GCMAMBOX documentation</title>
<style>
body {{ font-family: sans-serif; max-width: 900px; margin: 2em auto; padding: 0 1em; line-height: 1.5; }}
table {{ border-collapse: collapse; margin: 1em 0; }}
th, td {{ border: 1px solid #ccc; padding: 4px 10px; }}
th {{ background: #f0f0f0; }}
code, pre {{ background: #f5f5f5; }}
pre {{ padding: 10px; overflow-x: auto; }}
</style></head><body>{body}</body></html>"""

app.layout = html.Div([

#    html.Div([
#           html.Img(src='assets/logo.png', style={'height': '130px', 'width': '600px'}),],
#           style={'textAlign': 'center', 'marginBottom': '10px', 'padding': '10px'}),

html.Div([
    html.Div(style={'flex': '1'}),  # Spacer
    html.Img(src='assets/logo.png', style={'height': '130px', 'width': '419px'}),
    html.Div([
        html.A('📄 Doc !',
               href='/doc',
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
                html.Div("Set the parameters below manually, or select a predefined case to prefill them.",
                         style={'fontSize': '12px', 'marginBottom': '10px'}),

                # Run button at the top
                html.Button('Run Simulation', id='run-button', n_clicks=0,
                           style={'width': '100%', 'height': '40px', 'fontSize': '16px',
                                 'backgroundColor': '#4CAF50', 'color': 'white',
                                 'border': 'none', 'cursor': 'pointer', 'borderRadius': '5px', 'marginBottom': '10px'}),

                html.Div(id='status-message', style={'marginTop': '5px', 'marginBottom': '10px', 'fontWeight': 'bold', 'fontSize': '12px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Predefined cases: fill all parameters below from a namelist in assets/presets/ (values stay editable)
                html.H4("Predefined Case", style={'fontSize': '16px', 'marginTop': '10px'}),
                dcc.Dropdown(id='preset',
                             options=[{'label': 'Custom (keep current values)', 'value': 'custom'}] + list_presets(),
                             value='custom', clearable=False,
                             style={'fontSize': '12px', 'marginBottom': '5px'}),
                html.Div(id='preset-description', style={'fontSize': '11px', 'fontStyle': 'italic', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Time parameters
                html.H4("Time Parameters", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("Time step (s):", style={'fontSize': '12px'}),
                dcc.Input(id='mam_dt', type='number', value=600, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Number of steps:", style={'fontSize': '12px'}),
                dcc.Input(id='mam_nstep', type='number', value=100, style={'width': '100%', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Aerosol processes
                html.H4("Aerosol Processes", style={'fontSize': '16px', 'marginTop': '10px'}),
                dcc.Checklist(
                    id='process_controls',
                    options=[
                        {'label': ' Gas-Aerosol Exchange', 'value': 'gasaerexch'},
                        {'label': ' Rename', 'value': 'rename'},
                        {'label': ' New Nucleation', 'value': 'newnuc'},
                        {'label': ' Coagulation', 'value': 'coag'},
                    ],
                    value=['gasaerexch', 'rename', 'newnuc', 'coag'],
                    style={'fontSize': '12px'}
                ),

                html.Hr(style={'margin': '10px 0'}),

                # Meteorology
                html.H4("Meteorological Parameters", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("Temperature Initial (K):", style={'fontSize': '12px'}),
                dcc.Input(id='mtmin', type='number', value=285, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Temperature Final (K):", style={'fontSize': '12px'}),
                dcc.Input(id='mtmax', type='number', value=285, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Pressure (Pa):", style={'fontSize': '12px'}),
                dcc.Input(id='press', type='number', value=1.e5, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Relative Humidity Initial:", style={'fontSize': '12px'}),
                dcc.Input(id='mrhmin', type='number', value=0.4, step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Relative Humidity Final:", style={'fontSize': '12px'}),
                dcc.Input(id='mrhmax', type='number', value=0.6, step=0.01, style={'width': '100%', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Gas phase: initial concentrations + production rates
                # (SO2 and its oxidation are not used: mdo_gaschem = 0, qso2 = 0)
                html.H4("Gas Phase", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Div("Initial concentrations (ppb)", style={'fontSize': '13px', 'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Label("H2SO4:", style={'fontSize': '12px'}),
                dcc.Input(id='qh2so4', type='number', value=0., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("HNO3:", style={'fontSize': '12px'}),
                dcc.Input(id='qhno3', type='number', value=8., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("NH3:", style={'fontSize': '12px'}),
                dcc.Input(id='qnh3', type='number', value=10., style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("SOAG:", style={'fontSize': '12px'}),
                dcc.Input(id='qsoag', type='number', value=0.05, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("HCl:", style={'fontSize': '12px'}),
                dcc.Input(id='qhcl', type='number', value=0.02, style={'width': '100%', 'marginBottom': '10px'}),

                # Constant production rates of aerosol precursor gases (collapsible)
                html.Details([
                    html.Summary("Production rates (ppb/h)", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        html.Label("H2SO4:", style={'fontSize': '12px'}),
                        dcc.Input(id='h2so4_prod', type='number', value=0., step=0.001, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("HNO3:", style={'fontSize': '12px'}),
                        dcc.Input(id='hno3_prod', type='number', value=0., step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("NH3:", style={'fontSize': '12px'}),
                        dcc.Input(id='nh3_prod', type='number', value=0., step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("SOAG:", style={'fontSize': '12px'}),
                        dcc.Input(id='soag_prod', type='number', value=0., step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("HCl:", style={'fontSize': '12px'}),
                        dcc.Input(id='hcl_prod', type='number', value=0., step=0.01, style={'width': '100%', 'marginBottom': '5px'}),
                    ], style={'marginTop': '5px'})
                ], style={'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Aerosol initial conditions: total mass per mode + mass fractions
                html.H4("Aerosol Initial Conditions", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Div("Mode total mass (µg/m³)", style={'fontSize': '13px', 'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Label("Accum:", style={'fontSize': '12px'}),
                dcc.Input(id='masstot-0', type='number', value=10.0, min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Aitken:", style={'fontSize': '12px'}),
                dcc.Input(id='masstot-1', type='number', value=2.0,  min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Coarse:", style={'fontSize': '12px'}),
                dcc.Input(id='masstot-2', type='number', value=20.0, min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Prim C:", style={'fontSize': '12px'}),
                dcc.Input(id='masstot-3', type='number', value=1.5,  min=0, step=0.1, style={'width': '100%', 'marginBottom': '10px'}),

                # Mass fractions by mode
                html.Div("Mass fractions by mode", style={'fontSize': '13px', 'fontWeight': 'bold', 'marginBottom': '5px'}),

                # Mode 0
                html.Details([
                    html.Summary("Mode 0 - Accumulation", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 0, 0.30),
                        create_mass_fraction_input("POM", 0, 0.10),
                        create_mass_fraction_input("SOA", 0, 0.10),
                        create_mass_fraction_input("BC", 0, 0.05),
                        create_mass_fraction_input("DST", 0, 0.05),
                        create_mass_fraction_input("NCL", 0, 0.05),
                        create_mass_fraction_input("NO3", 0, 0.15),
                        create_mass_fraction_input("NH4", 0, 0.15),
                        create_mass_fraction_input("CO3", 0, 0.0),
                        create_mass_fraction_input("CA", 0, 0.0),
                        create_mass_fraction_input("CL", 0, 0.05),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 1
                html.Details([
                    html.Summary("Mode 1 - Aitken", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 1, 0.65),
                        create_mass_fraction_input("POM", 1, 0.0),
                        create_mass_fraction_input("SOA", 1, 0.20),
                        create_mass_fraction_input("BC", 1, 0.0),
                        create_mass_fraction_input("DST", 1, 0.0),
                        create_mass_fraction_input("NCL", 1, 0.0),
                        create_mass_fraction_input("NO3", 1, 0.0),
                        create_mass_fraction_input("NH4", 1, 0.15),
                        create_mass_fraction_input("CO3", 1, 0.0),
                        create_mass_fraction_input("CA", 1, 0.0),
                        create_mass_fraction_input("CL", 1, 0.0),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 2
                html.Details([
                    html.Summary("Mode 2 - Coarse", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        create_mass_fraction_input("SO4", 2, 0.0),
                        create_mass_fraction_input("POM", 2, 0.0),
                        create_mass_fraction_input("SOA", 2, 0.0),
                        create_mass_fraction_input("BC", 2, 0.0),
                        create_mass_fraction_input("DST", 2, 0.35),
                        create_mass_fraction_input("NCL", 2, 0.32),
                        create_mass_fraction_input("NO3", 2, 0.0),
                        create_mass_fraction_input("NH4", 2, 0.0),
                        create_mass_fraction_input("CO3", 2, 0.0),
                        create_mass_fraction_input("CA", 2, 0.0),
                        create_mass_fraction_input("CL", 2, 0.33),
                    ], style={'paddingLeft': '10px', 'paddingTop': '5px'})
                ], open=False, style={'marginBottom': '5px'}),

                # Mode 3
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

                html.Hr(style={'margin': '10px 0'}),

                # Accommodation coefficients (gas-aerosol exchange parameters, collapsible)
                html.Details([
                    html.Summary("Accommodation Coefficients", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        html.Label("H2SO4:", style={'fontSize': '12px'}),
                        dcc.Input(id='accom_coef_h2so4', type='number', value=0.650, min=0, max=1, step=0.001, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("HNO3:", style={'fontSize': '12px'}),
                        dcc.Input(id='accom_coef_hno3', type='number', value=0.100, min=0, max=1, step=0.001, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("HCl:", style={'fontSize': '12px'}),
                        dcc.Input(id='accom_coef_hcl', type='number', value=0.100, min=0, max=1, step=0.001, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("NH3:", style={'fontSize': '12px'}),
                        dcc.Input(id='accom_coef_nh3', type='number', value=0.092, min=0, max=1, step=0.001, style={'width': '100%', 'marginBottom': '5px'}),
                    ], style={'marginTop': '5px'})
                ], style={'marginBottom': '10px'}),


            ], style={'padding': '15px', 'overflowY': 'auto', 'minHeight': '100%'})
        ], style={'width': '15%', 'minHeight': '1650px', 'backgroundColor': '#f5f5f5', 'boxSizing': 'border-box'}),

        # Right panel - Visualization
        html.Div([
            dcc.Loading(
                id="loading",
                type="default",
                children=dcc.Graph(id='animation-plot', style={'height': '100%'})
            ),
        ], style={'width': '85%', 'height': '100%', 'padding': '10px', 'boxSizing': 'border-box', 'backgroundColor': '#e6f2ff'}),

    ], style={'display': 'flex', 'minHeight': '1650px'}),
], style={'fontFamily': 'Arial, sans-serif', 'margin': '0', 'padding': '0', 'backgroundColor': '#e6f2ff'})

# ============================================================================
# Callbacks - ONLY MODIFIED THE SIMULATION EXECUTION PART
# ============================================================================

@app.callback(
    [Output(field, 'value') for field in PRESET_FIELDS] + [Output('preset-description', 'children')],
    Input('preset', 'value'),
    prevent_initial_call=True
)
def load_preset(preset):
    """Fill the inputs from the selected preset file ('custom' leaves them unchanged)"""
    if preset == 'custom':
        return [dash.no_update] * len(PRESET_FIELDS) + ['']
    _, description, values = read_preset(os.path.join(PRESET_DIR, preset))
    return preset_to_inputs(values) + [description]

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
     State('qh2so4', 'value'),
     State('qhno3', 'value'),
     State('qnh3', 'value'),
     State('qsoag', 'value'),
     State('qhcl', 'value'),
     # Accommodation coefficients
     State('accom_coef_h2so4', 'value'),
     State('accom_coef_hno3', 'value'),
     State('accom_coef_hcl', 'value'),
     State('accom_coef_nh3', 'value'),
     # Gas production rates (ppb/h)
     State('h2so4_prod', 'value'),
     State('hno3_prod', 'value'),
     State('nh3_prod', 'value'),
     State('soag_prod', 'value'),
     State('hcl_prod', 'value'),
     # Total mass per mode (µg/m³)
     State('masstot-0', 'value'),
     State('masstot-1', 'value'),
     State('masstot-2', 'value'),
     State('masstot-3', 'value'),
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
                   qh2so4, qhno3, qnh3, qsoag, qhcl,
                   accom_coef_h2so4, accom_coef_hno3, accom_coef_hcl, accom_coef_nh3,
                   h2so4_prod, hno3_prod, nh3_prod, soag_prod, hcl_prod,
                   masstot0, masstot1, masstot2, masstot3,
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

    # For initial load, show prompt — never auto-load a stale mam_output.nc
    if n_clicks == 0:
        empty_fig = go.Figure().add_annotation(
            text="Check Parameters & Click 'Run Simulation' to start",
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
            'mdo_gaschem': 0,
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
            'masstot': [masstot0, masstot1, masstot2, masstot3],
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
            'qso2': 0.,
            'qh2so4': qh2so4,
            'qsoag': qsoag,
            'qhno3': qhno3,
            'qnh3': qnh3,
            'qhcl': qhcl,
            'accom_coef_h2so4': accom_coef_h2so4,
            'accom_coef_hno3': accom_coef_hno3,
            'accom_coef_hcl': accom_coef_hcl,
            'accom_coef_nh3': accom_coef_nh3,
            # empty input fields come back as None -> no source
            # h2so4_prod is in ppt/h in the namelist, ppb/h in the UI
            'h2so4_prod': (h2so4_prod or 0.) * 1000.,
            'hno3_prod': hno3_prod or 0.,
            'nh3_prod': nh3_prod or 0.,
            'soag_prod': soag_prod or 0.,
            'hcl_prod': hcl_prod or 0.,
        }

        # Copy executable to temp directory
        if os.path.exists(GCMAMBOX_EXE):
            shutil.copy(GCMAMBOX_EXE, temp_dir)
        else:
            return (go.Figure().add_annotation(
                text=f"gcmambox executable not found: {GCMAMBOX_EXE}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=20, color='red')
            ), "❌ Executable not found")

        # Copy optical input files to temp directory; gcmambox cannot run without all of them
        optical_files, missing_files = get_optical_input_files(ASSETS_DIR)
        if missing_files:
            return (go.Figure().add_annotation(
                text=f"Optical files missing in {ASSETS_DIR}:<br>" + "<br>".join(missing_files),
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16, color='red')
            ), "❌ Optical input files missing")
        for fpath in optical_files:
            shutil.copy(fpath, temp_dir)


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
            fig = create_animated_figure(session_output, dt=mam_dt, press=press)

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
    import argparse
    parser = argparse.ArgumentParser(description='MAM box model dashboard')
    parser.add_argument('-p', '--port', type=int, default=8050, help='port to serve on (default: 8050)')
    args = parser.parse_args()
    app.run(debug=True, host='0.0.0.0', port=args.port)

