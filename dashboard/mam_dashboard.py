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
GCMAMBOX_EXE = os.path.join(SCRIPT_DIR, "..", "modbin", 'gcmambox')
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

def coltr(colo, t):
    """Named color -> 'rgba(r, g, b, t)' with 0-255 channels, as plotly expects."""
    r, g, b = (round(255 * v) for v in mpc.to_rgb(colo))
    return f'rgba({r}, {g}, {b}, {t})'

# Transparency of the dm/dlogD fills (only water is transparent). It must be in
# the fill color: with fill='tonexty' plotly draws each fill in the previous
# trace's layer, so the trace 'opacity' does not apply to its own fill.
tr = 1.0
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
'wat_aer': coltr('lightblue', 0.15),

    }

modname = {0: 'Accumulation', 1: 'Aitken', 2: 'Coarse', 3: 'Primary Carbon'}

sigmag = np.array([1.800, 1.600, 1.800, 1.600])

# Mode colours, shared by the number panel and the mode-property panels; chosen apart
# from the species colours. modfill = same colour, semi-transparent, for overlaid fills
modcolor = {0: '#2a9d8f',    # Accumulation: blue-green (teal)
            1: '#a78bfa',    # Aitken: light violet
            2: '#e9894a',    # Coarse: soft orange
            3: '#8b5a3c'}    # Primary Carbon: brown
MODE_FILL_ALPHA = 0.45
modfill = {m: 'rgba' + str(tuple(int(c[i:i+2], 16) for i in (1, 3, 5)) + (MODE_FILL_ALPHA,))
           for m, c in modcolor.items()}
MODE_ORDER = [1, 3, 0, 2]   # Aitken, Primary Carbon, Accumulation, Coarse (by size)
MODE_SHORT = {0: 'Accum.', 1: 'Aitken', 2: 'Coarse', 3: 'PC'}
PROP_VARS = ['aod_mode', 'ssa_mode', 'asm_mode', 'aod_lw_mode', 'hygro_aer', 'pH_aer']
PROP_TITLES = ['AOD 550 nm', 'Single scattering albedo', 'Asymmetry parameter g (light →)',
               'LW absorption optical depth, 9.3–10.2 µm', 'Hygroscopicity κ', 'pH']
PH_RANGE = [-1, 8]
MAX_FRAMES = 150      # animation frames per run (steps are thinned beyond that)
MAX_STEPS = 2000      # longest run accepted (server load)
MAX_SUBSTEPS = 10     # MOSAIC/amicphys substeps per time step (mam_nsubstep)
MIN_MODE_MASS = 1.e-3   # µg/m³
# Phase-function lobes of the asymmetry panel: half-size in x (mode units) and y (g units)
LOBE_THETA = np.linspace(0., 2. * np.pi, 61)
LOBE_RX, LOBE_RY = 0.42, 0.14
G_EMPTY = 0.5
# Colour scales of the circles (SSA, hygroscopicity, pH), with the values mapped to
# both ends. 'colorscale' can also be any plotly named scale, e.g. 'Viridis', 'Plasma',
# 'Teal', 'ice', 'gray' (add '_r' to reverse): see plotly.express.colors.sequential
CIRCLE_SCALES = {
    'ssa_mode':  dict(colorscale=[[0, 'black'], [1, '#f8f8f8']], cmin=0.5, cmax=1.0),
    # yellow (weakly hygroscopic) -> pure blue, through near-white rather than the
    # grey/green of a direct yellow-blue blend; still light at kappa ~ 0.4
    'hygro_aer': dict(colorscale=[[0, '#ffe066'], [0.3, '#f5f8ff'], [0.6, '#7fb2ff'],
                                  [1, '#0033ff']], cmin=0.0, cmax=1.0),
    # magenta (very acidic) -> violet -> green -> yellow at neutral pH 7; no blue step
    # (kept for hygroscopicity), and a short violet-green transition to limit the grey
    'pH_aer':    dict(colorscale=[[0, '#e040fb'], [0.36, '#9c6cff'], [0.4, '#66d17a'],
                                  [1, '#ffee58']], cmin=0.0, cmax=7.0),
}

# ============================================================================
# Utility functions - ONLY MODIFIED write_namelist and run_fortran_model
# ============================================================================

def write_namelist(params, namelist_path='namelist'):
    """Write namelist file from parameters dictionary"""
    # written only when > 1, so the namelist still works with binaries built before mam_nsubstep
    nsubstep = int(params.get('mam_nsubstep', 1) or 1)
    nsubstep_line = f"\nmam_nsubstep   = {nsubstep}," if nsubstep > 1 else ""
    namelist_content = f"""&time_input
mam_dt         = {params['mam_dt']},
mam_nstep      = {params['mam_nstep']},{nsubstep_line}
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
PRESET_SCALARS = ['mam_dt', 'mam_nstep', 'mam_nsubstep', 'mtmin', 'mtmax', 'press', 'mrhmin', 'mrhmax',
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
            out[key] = int(values[key]) if key in ('mam_dt', 'mam_nstep', 'mam_nsubstep') else values[key]
    if 'mam_dt' in values:
        out.setdefault('mam_nsubstep', 1)     # presets without it: back to standard MAM
    for m, v in enumerate(values.get('masstot', [])):
        out[f'masstot-{m}'] = v
    for sp in MF_SPECIES:
        for m, v in enumerate(values.get(f'mf{sp}', [])):
            out[f'mf{sp}-{m}'] = v
    return [out.get(field, dash.no_update) for field in PRESET_FIELDS]

def _round_sig(v, sig=4):
    """Round numbers (array or list) to `sig` significant digits, NaN -> None.
    Vectorized: round(a * 10**d) / 10**d gives floats with a short repr in JSON"""
    if isinstance(v, (list, tuple)) and any(isinstance(x, str) for x in v):
        return v                                   # category labels
    try:
        a = np.asarray(v, dtype=float)
    except (TypeError, ValueError):
        return v
    with np.errstate(divide='ignore', invalid='ignore'):
        dec = sig - 1 - np.floor(np.log10(np.abs(a)))
    dec = np.where(np.isfinite(dec), dec, 0)
    scale = 10.0 ** np.abs(dec)
    r = np.where(dec >= 0, np.round(a * scale) / scale, np.round(a / scale) * scale)
    out = r.tolist()
    if np.isnan(a).any():                          # None (JSON null) for NaN
        out = np.where(np.isnan(a), None, r.astype(object)).tolist()
    return out

class _StepView:
    """Fast in-memory stand-in for ds.sel(nsteps=s) in create_animated_figure:
    df[var].sel(mode=m), .values, float(), df.mode.values, df[var].sum()"""
    def __init__(self, arrays, modes, s):
        self._arrays, self._s = arrays, s
        self.mode = _VarView(modes, {})
        self._modes = {int(m): i for i, m in enumerate(modes)}

    def __getitem__(self, var):
        return _VarView(self._arrays[var][self._s], self._modes)

class _VarView:
    def __init__(self, a, modes):
        self.values, self._modes = a, modes

    def sel(self, mode):
        return _VarView(self.values[self._modes[int(mode)]], self._modes)

    def sum(self):
        return _VarView(self.values.sum(), self._modes)

    def __float__(self):
        return float(self.values)

    def __mul__(self, other):
        return self.values * other
    __rmul__ = __mul__

def _trace(kind, **kw):
    """Plotly trace as a plain dict. Much faster than go.Scatter/go.Bar for the
    hundreds of animation frames, which plotly would otherwise validate one by one.
    Expands the marker_/line_ shortcuts and rounds the data arrays to 4 significant
    digits to keep the page small."""
    out = {'type': kind}
    for k, v in kw.items():
        if k in ('x', 'y', 'base', 'customdata'):
            v = _round_sig(v)
        prefix, _, sub = k.partition('_')
        if sub and prefix in ('marker', 'line'):
            out.setdefault(prefix, {})[sub] = v
        else:
            out[k] = v
    return out

def create_animated_figure(nc_file='mam_output.nc', dt=1200, press=1.e5):
    """Create the animated plotly figure from NetCDF output"""
    try:
        # load into memory once: reading values lazily from the file is slow
        with xr.open_dataset(nc_file) as f:
            ds = f.load()
        print('open', nc_file)
    except FileNotFoundError:
        return go.Figure().add_annotation(
            text="No simulation data available. Run simulation first.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )
    maxyaod = ds['aod_mode'].max(dim='nsteps').max(dim='mode').values * 1.4

    # Animate at most MAX_FRAMES steps (every k-th step, always including the last);
    # the gas time series still use every step
    nsteps = len(ds.nsteps)
    stride = -(-nsteps // MAX_FRAMES)
    frame_steps = [int(s) for s in ds.nsteps.values[::stride]]
    if frame_steps[-1] != int(ds.nsteps.values[-1]):
        frame_steps.append(int(ds.nsteps.values[-1]))
    t_end = nsteps * dt / 3600.   # run length in hours

    # per-variable numpy arrays (nsteps first): much faster than xarray .sel per value
    arrays = {v: ds[v].transpose('nsteps', ...).values
              for v in ds.data_vars if 'nsteps' in ds[v].dims}
    modes = ds['mode'].values

    frames = []
    max_dmdlogd = 0.   # mass panel y range: max over all frames, not the last one
    for s in frame_steps:
        df = _StepView(arrays, modes, s)
        data_list = []
        ntrac1 = 0
        # Gas concentrations
        for g in list(gascolor.keys()):
            conc = ds[g].sel(nsteps=slice(0, int(s+1)))
            trace = _trace('scatter',
                x=conc.nsteps.values * dt / 3600.,
                y=conc.values,
                # dot on the current value only, so the curve is visible from step 0
                mode='lines+markers',
                marker=dict(size=[0] * (conc.size - 1) + [7]),
                name=f'{g} ({float(df[g].values):.3g})',
                line_color=gascolor[g],
                legend='legend'
            )
            data_list.append(trace)
            ntrac1 = ntrac1 +1

        # RH and TEMP gauges — top-right of subplot 1 (gas concentrations)
        rh_value   = float(df['relhum'] * 100)
        temp_value = float(df['temp'])

        x_rh_pos   = t_end * 0.78
        x_temp_pos = t_end * 0.92
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

        trace_rh_bg = _trace('scatter',
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

        trace_rh_text = _trace('scatter',
            x=[x_rh_pos], y=[y_rh_pos],
            mode='text',
            text=[f'<b>RH</b><br>{rh_value:.0f}<br><b>%</b>'],
            textfont=dict(size=12, color=text_color_rh, family='Arial Black'),
            textposition='middle center',
            showlegend=False, hoverinfo='skip', name=''
        )
        data_list.append(trace_rh_text)
        ntrac1 = ntrac1 + 1

        trace_temp_bg = _trace('scatter',
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

        trace_temp_text = _trace('scatter',
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
            trace_name = f'{c} ({sp_total:.3g})' if sp_total > 1.e-3 else c

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

            # water: lighter and hatched, so it reads as water and not as dry mass
            water = c == 'wat_aer'
            trace = _trace('scatter',
                x=x_vol, y=ycum, fill=fil,
                mode='none', fillcolor=specolor[c],
                name=trace_name,
                legend='legend2',
                customdata=cdata, hovertemplate=htemplate,
                **({'fillpattern': {'shape': '/', 'fgcolor': 'steelblue', 'size': 8, 'solidity': 0.3}}
                   if water else {})
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
        trace = _trace('scatter',
            x=x_vol, y=ycum, mode='lines',
            name=f'Total Mass ({grand_total:.2f})', line_color='grey',
            legend='legend2',
            customdata=total_cdata, hovertemplate=total_htemplate
        )
        data_list.append(trace)
        ntrac2 = ntrac2+1
        max_dmdlogd = max(max_dmdlogd, float(np.max(ycum)))

# Number distribution, converted from #/kg-air to #/cm3 with the air density at this step
        kg_to_cm3 = press / (287.05 * float(df['temp'])) * 1.e-6
        ymodcum = 0
        ntrac3 = ntrac2
        for m in df.mode.values:
            Dg = df['Dgn_mode'].sel(mode=m).values * 1.E6
            sigma_g = sigmag[m]
            num = df['num_aer'].sel(mode=m).values * kg_to_cm3
            x_num, y_num = lognormal_distribution(Dg, sigma_g, num)
            trace = _trace('scatter',
                x=x_num, y=y_num, fill='tozeroy',
                mode='none', name=f'{modname[int(m)]} ({float(num):.3g})',
                line_color=modcolor[int(m)], fillcolor=modfill[int(m)], legend='legend3'
            )
            data_list.append(trace)
            ymodcum = ymodcum + y_num
            ntrac3 = ntrac3+1
        trace = _trace('scatter',
            x=x_num, y=ymodcum, mode='lines',
            name=f"Total number ({float(df['num_aer'].sum().values) * kg_to_cm3:.3g})", line_color='black',
            legend='legend3'
        )
        data_list.append(trace)
        ntrac3 = ntrac3+1

#  Mode properties: 2 x 3 small bar charts (AOD by species, SSA, g, LW abs, kappa, pH)
        ntrac4 = ntrac3
        p4_bounds = []          # trace index after each mini panel, in PROP_PANELS order
        x_positions = [MODE_SHORT[m] for m in MODE_ORDER]
        # modes with negligible mass (e.g. 1e-13 µg/m³ left by rename) are shown empty
        mode_has_mass = [mode_total_mass.get(m, 0.0) > MIN_MODE_MASS for m in MODE_ORDER]

        def mode_values(var):
            if var not in ds.variables:
                return [np.nan] * len(MODE_ORDER)
            return [float(df[var].sel(mode=m).values) for m in MODE_ORDER]

        props = {
            'aod_mode':    mode_values('aod_mode'),
            'ssa_mode':    mode_values('ssa_mode'),
            'asm_mode':    mode_values('asm_mode'),
            'aod_lw_mode': mode_values('aod_lw_mode'),
            'hygro_aer':   mode_values('hygro_aer'),
            'pH_aer':      mode_values('pH_aer'),
        }
        # Same hover box on every bar: all the properties of the mode
        cdata = [[props[v][i] for v in PROP_VARS] for i in range(len(MODE_ORDER))]
        htemplate = ('<b>%{x}</b><br>'
                     'AOD 550: %{customdata[0]:.3g}<br>SSA: %{customdata[1]:.3f}<br>'
                     'g: %{customdata[2]:.3f}<br>LW abs.: %{customdata[3]:.3g}<br>'
                     'κ: %{customdata[4]:.3f}<br>pH: %{customdata[5]:.2f}<extra></extra>')

        # (a) AOD at 550 nm, stacked by species
        aod_species = ['aod_sulfate', 'aod_bc', 'aod_pom', 'aod_mom', 'aod_soa',
                       'aod_dust', 'aod_seasalt', 'aod_no3', 'aod_nh4']
        species_colors = {
            'aod_sulfate': 'green', 'aod_bc': 'black', 'aod_pom': 'magenta',
            'aod_mom': 'brown', 'aod_soa': 'indigo', 'aod_dust': 'yellow',
            'aod_seasalt': 'grey', 'aod_no3': 'blue', 'aod_nh4': 'orange',
        }
        for species in aod_species:
            if species in ds.variables:
                data_list.append(_trace('bar',
                    x=x_positions, y=mode_values(species),
                    name=species.replace('aod_', ''),
                    marker_color=species_colors[species], opacity=0.75,
                    legend='legend4',
                    customdata=cdata, hovertemplate=htemplate,
                ))
                ntrac4 += 1
        # total AOD written on top of each stack: an invisible zero-height bar stacked
        # last, so the label moves with the bars during the animation
        data_list.append(_trace('bar',
            x=x_positions, y=[1.e-12] * len(MODE_ORDER), marker_opacity=0,
            text=[f'{v:.3g}' if ok else '' for v, ok in zip(props['aod_mode'], mode_has_mass)],
            textposition='outside', cliponaxis=False,
            showlegend=False, hoverinfo='skip',
        ))
        ntrac4 += 1
        p4_bounds.append(ntrac4)

        def circles(var, fmt, empty_label, y_empty):
            # One circle per mode, coloured by value. Every mode always keeps its point
            # (fixed ids, invisible when empty) so that plotly animates each circle up
            # and down instead of shuffling the points when a mode appears/disappears.
            values = masked(var)
            scale = CIRCLE_SCALES[var]
            empty = [v is None or np.isnan(v) for v in values]
            return _trace('scatter',
                x=x_positions, y=[y_empty if e else v for v, e in zip(values, empty)],
                ids=x_positions, mode='markers+text',
                marker=dict(size=24, color=[scale['cmin'] if e else v for v, e in zip(values, empty)],
                            **scale,
                            opacity=[0 if e else 0.75 for e in empty],
                            line=dict(color='black', width=1)),
                text=[empty_label if e else format(v, fmt) for v, e in zip(values, empty)],
                textposition=['middle center' if e else 'middle right' for e in empty],
                textfont=dict(size=12, color=['grey' if e else '#444' for e in empty]),
                showlegend=False, customdata=cdata, hovertemplate=htemplate,
                cliponaxis=False,
            )

        def masked(var):
            return [v if ok else None for v, ok in zip(props[var], mode_has_mass)]

        # (b) SSA: circles, dark (absorbing) -> light (scattering)
        data_list.append(circles('ssa_mode', '.3f', '', 0.5))
        ntrac4 += 1
        p4_bounds.append(ntrac4)

        # (c) asymmetry parameter: Henyey-Greenstein phase-function lobe of each mode,
        # placed at height g, light coming from the left: round for g = 0 (isotropic),
        # stretched to the right for forward scattering. Radius ~ (P/Pmax)^(1/3) so
        # that the lobe stays visible for large g. Fixed number of points per lobe and
        # fixed ids, so that the lobes animate smoothly; empty modes collapse to a point.
        g_vals = [None if v is None or np.isnan(v) else v for v in masked('asm_mode')]
        for i, g in enumerate(g_vals):
            if g is None:
                xs = np.full(LOBE_THETA.size, float(i))
                ys = np.full(LOBE_THETA.size, G_EMPTY)
            else:
                p = (1 - g**2) / (1 + g**2 - 2 * g * np.cos(LOBE_THETA))**1.5
                r = (p / p.max())**(1. / 3.)
                xs = i + LOBE_RX * r * np.cos(LOBE_THETA)
                ys = g + LOBE_RY * r * np.sin(LOBE_THETA)
            data_list.append(_trace('scatter',
                x=xs, y=ys, mode='lines', fill='toself',
                line=dict(color=modcolor[MODE_ORDER[i]], width=1),
                fillcolor=modfill[MODE_ORDER[i]],
                showlegend=False, hoverinfo='skip',
            ))
        # the particle (dot, with hover) and the value written above the lobe
        data_list.append(_trace('scatter',
            x=list(range(len(MODE_ORDER))), y=[G_EMPTY if g is None else g for g in g_vals],
            ids=x_positions, mode='markers',
            marker=dict(size=6, color='black', opacity=[0 if g is None else 1 for g in g_vals]),
            showlegend=False, customdata=cdata, hovertemplate=htemplate,
        ))
        data_list.append(_trace('scatter',
            x=list(range(len(MODE_ORDER))),
            y=[G_EMPTY if g is None else g + LOBE_RY for g in g_vals],
            ids=x_positions, mode='text',
            text=['' if g is None else f'{g:.3f}' for g in g_vals],
            textposition='top center', textfont=dict(size=12, color='#444'),
            showlegend=False, hoverinfo='skip', cliponaxis=False,
        ))
        ntrac4 += len(MODE_ORDER) + 2
        p4_bounds.append(ntrac4)

        # (d) LW absorption: one bar per mode, mode colours
        yv = masked('aod_lw_mode')
        data_list.append(_trace('bar',
            x=x_positions, y=yv,
            marker_color=[modcolor[m] for m in MODE_ORDER], opacity=0.8,
            text=[format(v, '.2g') if v is not None else '' for v in yv],
            textposition='outside', cliponaxis=False,
            showlegend=False, customdata=cdata, hovertemplate=htemplate,
        ))
        ntrac4 += 1
        p4_bounds.append(ntrac4)

        # (e) hygroscopicity: circles
        data_list.append(circles('hygro_aer', '.2f', '', 0.6))
        ntrac4 += 1
        p4_bounds.append(ntrac4)

        # (f) pH: circles; "no water" written mid-panel where pH is undefined
        data_list.append(circles('pH_aer', '.2f', 'no water', PH_RANGE[1] / 2))
        ntrac4 += 1
        p4_bounds.append(ntrac4)

        frames.append(dict(name='step%s' % s, data=data_list))

    # Create figure: 3 full-width panels, then the mode properties as 2 x 3 small panels
    full = [{'colspan': 3}, None, None]
    fig = go.Figure().set_subplots(
        5, 3, specs=[full, full, full, [{}, {}, {}], [{}, {}, {}]],
        vertical_spacing=0.05, horizontal_spacing=0.07,
        row_heights=[0.20, 0.24, 0.20, 0.18, 0.18],
        subplot_titles=('<b>Gas concentrations</b>',
                        '<b>Aerosol mass size distribution</b>',
                        '<b>Aerosol number size distribution</b>')
                       + tuple(f'<b>{t}</b>' for t in PROP_TITLES))

    fig.add_traces(frames[0]['data'][0:ntrac1], 1, 1)
    fig.add_traces(frames[0]['data'][ntrac1:ntrac2], 2, 1)
    fig.add_traces(frames[0]['data'][ntrac2:ntrac3], 3, 1)
    start = ntrac3
    for k, end in enumerate(p4_bounds):
        fig.add_traces(frames[0]['data'][start:end], 4 + k // 3, 1 + k % 3)
        start = end

    # Update axes
    fig.update_xaxes(range=[-0.01 * t_end, t_end], title='time (h)', row=1, col=1)
    fig.update_yaxes(type='log', range=[-2., np.log10(70.)], title='gas conc. (ppb)', row=1, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=2, col=1)
    fig.update_yaxes(range=[0, 1.05 * max_dmdlogd or 1.], title='aer. dm/dlogD (µg.m⁻³)', row=2, col=1)

    fig.update_xaxes(type="log", range=[-3, 2], title='D (microm)', row=3, col=1)
    fig.update_yaxes(type='log', range=[-1, 6], title='aer. dN/dlogD (#.cm⁻³)', row=3, col=1)

    # Mode properties: fixed ranges so that the animation does not rescale
    maxlw = float(ds['aod_lw_mode'].max()) * 1.3 if 'aod_lw_mode' in ds.variables else 1.0
    prop_ranges = [[0, maxyaod], [0, 1.1], [-0.05, 1.25], [0, maxlw or 1.0], [-0.1, 1.4], PH_RANGE]
    for k, yr in enumerate(prop_ranges):
        fig.update_yaxes(range=yr, row=4 + k // 3, col=1 + k % 3)
        fig.update_xaxes(tickfont=dict(size=11), row=4 + k // 3, col=1 + k % 3)
    fig.update_yaxes(title='wet, 1 km layer', row=4, col=1)
    # asymmetry panel uses a numeric x axis (lobes drawn around each mode position)
    fig.update_xaxes(tickvals=list(range(len(MODE_ORDER))),
                     ticktext=[MODE_SHORT[m] for m in MODE_ORDER],
                     range=[-0.6, len(MODE_ORDER) - 0.4], showgrid=False, zeroline=False,
                     row=4, col=3)


# Animation controls
    fr_duration = 500
    sliders = [{
        "pad": {"b": 10, "t": 10},
        "len": 0.9,
        "x": 0.1,
        "y": 1.025,  # above the first subplot title
        "yanchor": "bottom",
        "currentvalue":{"prefix": "Time (h) : "},
        "steps": [
            {
                "args": [[f['name']], {
                    "frame": {"duration": fr_duration},
                    "mode": "immediate",
                    "fromcurrent": True,
                    "transition": {"duration": fr_duration, "easing": "linear"},
                }],
#                "label": f"tstep{k+1}",add time label in mn
                "label": f'{s * dt / 3600.:.2f}',
                "method": "animate",
            }
            for s, f in zip(frame_steps, frames)
        ],
    }]

    fig.update_layout(
        sliders=sliders,
        barmode='stack',  # Enable stacking for AOD bars
        height=1800,
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
            traceorder='normal',
            title=dict(text='<b>AOD by species</b>', font=dict(size=14))
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

    # Frames are added as plain dicts (no plotly validation, see _trace); dcc.Graph
    # accepts the figure as a dict
    fig_dict = fig.to_dict()
    fig_dict['frames'] = frames
    return fig_dict

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

app = dash.Dash(__name__, title="GCMAMBOX")

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
        script = ""
    except ImportError:
        # No python markdown package (e.g. on a server): render in the browser
        # with marked.js; the plain text stays visible if the CDN is unreachable
        body = f'<pre id="doc-src">{htmllib.escape(text)}</pre>'
        script = """<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<script>
if (window.marked) {
  const src = document.getElementById('doc-src');
  const div = document.createElement('div');
  div.innerHTML = marked.parse(src.textContent);
  src.replaceWith(div);
}
</script>"""
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>GCMAMBOX documentation</title>
<style>
body {{ font-family: sans-serif; max-width: 900px; margin: 2em auto; padding: 0 1em; line-height: 1.5; }}
table {{ border-collapse: collapse; margin: 1em 0; }}
th, td {{ border: 1px solid #ccc; padding: 4px 10px; }}
th {{ background: #f0f0f0; }}
code, pre {{ background: #f5f5f5; }}
pre {{ padding: 10px; overflow-x: auto; }}
</style></head><body>{body}{script}</body></html>"""

def serve_layout():
    """Page layout, rebuilt at each page load so the cases menu sees new or removed preset files"""
    return html.Div([

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
                             options=list_presets(), value=None, placeholder='Choose a case...',
                             clearable=False,
                             style={'fontSize': '12px', 'marginBottom': '5px'}),
                html.Div(id='preset-description', style={'fontSize': '11px', 'fontStyle': 'italic', 'marginBottom': '10px'}),

                html.Hr(style={'margin': '10px 0'}),

                # Time parameters
                html.H4("Time Parameters", style={'fontSize': '16px', 'marginTop': '10px'}),
                html.Label("Time step (s):", style={'fontSize': '12px'}),
                dcc.Input(id='mam_dt', type='number', value=300, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("Number of steps:", style={'fontSize': '12px'}),
                dcc.Input(id='mam_nstep', type='number', value=200, style={'width': '100%', 'marginBottom': '5px'}),
                html.Label("MOSAIC substeps:", style={'fontSize': '12px'}),
                dcc.Input(id='mam_nsubstep', type='number', value=1, min=1, max=MAX_SUBSTEPS, step=1,
                          style={'width': '100%', 'marginBottom': '10px'}),

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
                html.Details([
                    html.Summary("Initial concentrations (ppb)", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
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
                    ], style={'marginTop': '5px'})
                ], style={'marginBottom': '10px'}),

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
                html.Details([
                    html.Summary("Mode total mass (µg/m³)", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
                    html.Div([
                        html.Label("Accum:", style={'fontSize': '12px'}),
                        dcc.Input(id='masstot-0', type='number', value=10.0, min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("Aitken:", style={'fontSize': '12px'}),
                        dcc.Input(id='masstot-1', type='number', value=2.0,  min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("Coarse:", style={'fontSize': '12px'}),
                        dcc.Input(id='masstot-2', type='number', value=20.0, min=0, step=0.1, style={'width': '100%', 'marginBottom': '5px'}),
                        html.Label("Prim C:", style={'fontSize': '12px'}),
                        dcc.Input(id='masstot-3', type='number', value=1.5,  min=0, step=0.1, style={'width': '100%', 'marginBottom': '10px'}),
                    ], style={'marginTop': '5px'})
                ], style={'marginBottom': '10px'}),

                # Mass fractions by mode
                html.Div("Mass fractions by mode", style={'fontSize': '13px', 'fontWeight': 'bold', 'marginTop': '20px', 'marginBottom': '5px'}),

                # Mode 0
                html.Details([
                    html.Summary("Mode 1 - Accumulation", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
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
                    html.Summary("Mode 2 - Aitken", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
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
                    html.Summary("Mode 3 - Coarse", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
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
                    html.Summary("Mode 4 - Primary Carbon", style={'fontSize': '13px', 'fontWeight': 'bold', 'cursor': 'pointer'}),
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
        ], style={'width': '15%', 'backgroundColor': '#f5f5f5', 'boxSizing': 'border-box'}),

        # Right panel - Visualization
        html.Div([
            dcc.Loading(
                id="loading",
                type="default",
                children=dcc.Graph(id='animation-plot', style={'height': '100%'})
            ),
        ], style={'width': '85%', 'height': '100%', 'padding': '10px', 'boxSizing': 'border-box', 'backgroundColor': '#e6f2ff'}),

    ], style={'display': 'flex'}),
], style={'fontFamily': 'Arial, sans-serif', 'margin': '0', 'padding': '0', 'backgroundColor': '#e6f2ff'})

app.layout = serve_layout

# ============================================================================
# Callbacks - ONLY MODIFIED THE SIMULATION EXECUTION PART
# ============================================================================

@app.callback(
    [Output(field, 'value') for field in PRESET_FIELDS] + [Output('preset-description', 'children')],
    Input('preset', 'value'),
    prevent_initial_call=True
)
def load_preset(preset):
    """Fill the inputs from the selected preset file (none selected: inputs unchanged)"""
    if not preset:
        return [dash.no_update] * len(PRESET_FIELDS) + ['']
    _, description, values = read_preset(os.path.join(PRESET_DIR, preset))
    return preset_to_inputs(values) + [description]

@app.callback(
    [Output('animation-plot', 'figure'),
     Output('status-message', 'children')],
    [Input('run-button', 'n_clicks')],
    [State('mam_dt', 'value'),
     State('mam_nstep', 'value'),
     State('mam_nsubstep', 'value'),
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
def run_simulation(n_clicks, mam_dt, mam_nstep, mam_nsubstep, processes, mtmin,mtmax, press, mrhmin,mrhmax,
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

    # Keep runs short on a shared server
    if not mam_nstep or mam_nstep > MAX_STEPS:
        message = f"Number of steps must be between 1 and {MAX_STEPS} (got {mam_nstep})"
        error_fig = go.Figure().add_annotation(
            text=message, xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False, font=dict(size=20, color='red'))
        return error_fig, f"❌ {message}"
    if not mam_nsubstep or not 1 <= mam_nsubstep <= MAX_SUBSTEPS:
        message = f"MOSAIC substeps must be between 1 and {MAX_SUBSTEPS} (got {mam_nsubstep})"
        error_fig = go.Figure().add_annotation(
            text=message, xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False, font=dict(size=20, color='red'))
        return error_fig, f"❌ {message}"

    # Create unique temporary directory for this simulation
    session_id = str(uuid.uuid4())
    temp_dir = tempfile.mkdtemp(prefix=f'mam_sim_{session_id}_')

    try:
        # Build parameters dictionary
        params = {
            'mam_dt': mam_dt,
            'mam_nstep': mam_nstep,
            'mam_nsubstep': int(mam_nsubstep),
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
            print('mam_dt',mam_dt)
            # Create figure from session-specific results
            fig = create_animated_figure(output_nc, dt=mam_dt, press=press)

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
    app.run(debug=True, dev_tools_ui=False, host='0.0.0.0', port=args.port)

