# GCMAMBOX — MAM aerosol box model

GCMAMBOX runs the **Modal Aerosol Module (MAM4)** in a single air parcel (box
model). It follows how an aerosol population and its precursor gases evolve
over time under prescribed temperature, humidity and pressure, and computes
other key modal aerosol properties such as optical properties,
hygroscopicity and pH. The dashboard lets you set up a case, run
it, and watch the evolution as an animation.

**Contact:** Fabien Solmon, Atmospheric Chemistry group, IAC, ETH Zürich —
fsolmon -at- env.ethz.ch. Code: <https://github.com/fsolmon/GCMAMBOX>.

### Origin of the box model

GCMAMBOX is built on a MAM4-MOSAIC library: the Modal Aerosol Model (MAM;
Liu et al., 2012) interacting with the Model for Simulating Aerosol
Interactions and Chemistry (MOSAIC; Zaveri et al., 2008). This library is
derived from the E3SM-MOSAIC code version of Wu et al. (2022), and it is the
same library that is coupled with the GEOS-Chem chemical transport model.

---

## 1. The model

### Aerosol representation

The aerosol is described by **4 log-normal modes**. Each mode has a number
concentration, a mass per chemical species, and aerosol water. At start-up
each mode gets its default dry diameter; afterwards the diameter evolves with
the mode's mass and number.

| # | Mode | Default dry diameter | σ<sub>g</sub> |
|---|------|---------------------:|-----:|
| 1 | Accumulation | 0.110 µm | 1.8 |
| 2 | Aitken | 0.026 µm | 1.6 |
| 3 | Coarse | 2.0 µm | 1.8 |
| 4 | Primary carbon | 0.050 µm | 1.6 |

Species: sulfate (SO4), primary organic matter (POM), secondary organic
aerosol (SOA), black carbon (BC), dust (DST non CaCO3), nitrate (NO3),
ammonium (NH4), carbonate (CO3), calcium (Ca), sodium(Na), chloride (Cl). Aerosol water mass is calculated in function of RH and modal properties. Not every mode
carries every species; a mass fraction given for a species a mode does not
carry is ignored.

Gases: H2SO4, HNO3, NH3, HCl and SOAG (the gas-phase SOA precursor).


### Processes

| Process | What it does |
|---|---
| **Gas–aerosol exchange** | Condensation/evaporation of H2SO4, HNO3, HCl, NH3 and SOAG on each mode. Inorganic partitioning and aerosol thermodynamics are computed with **MOSAIC** (it includes nitrate, ammonium and chloride, plus aerosol pH) |
| **Rename** | Moves particles that have grown past the Aitken size range into the accumulation mode |
| **New particle formation** | Nucleation of new particles from H2SO4 (and NH3) into the Aitken mode |
| **Coagulation** | Collisions between particles, within and between modes |
| **Primary-carbon aging** | Moves primary-carbon particles (fresh BC/POM) to the accumulation mode once they are coated by a few molecular layers of condensed material or by coagulation; part of gas–aerosol exchange and coagulation |
| **Aitken ↔ accumulation size adjustment** | Always on (the rename switch does not turn it off). Each mode is a lognormal of fixed width with a size range; when the median size of the Aitken or accumulation mode leaves its range, part of the particles is moved to the other mode. The Aitken mode can therefore receive particles that were not nucleated (see the Black carbon aging case) |
| **Gas sources** | No explicit gas-phase chemistry: the precursor gases (H2SO4, HNO3, NH3, HCl, SOAG) come only from their initial concentrations andi/or from constant production rates set by the user (ppb/h), added at every step |

The four aerosol processes (gas–aerosol exchange, rename, nucleation,
coagulation) can each be switched on or off.

Water uptake (hygroscopic growth) is always computed from the ambient RH and
the hygroscopicity of each mode's composition. There are no clouds in the box for now, although there is a specific cloud path in MAM ( maybe added to dashboard in the future) .

### Optical properties

At every step the model computes, from the wet size and composition of each mode:

- **Shortwave**: aerosol optical depth (AOD), single scattering albedo and
  asymmetry parameter, per mode and per species, on the RRTMG bands
  (the dashboard shows AOD at 550 nm).
- **Longwave**: absorption optical depth on the RRTMG longwave bands
  (the dashboard shows the 9.3–10.2 µm band, in the atmospheric window).

The parameterisation is the CESM (CAM) modal optics of Ghan and Zaveri
(2007): the refractive index of each mode is the volume average of its
components, water included, and the optical properties come from Chebyshev
fits of Mie calculations as a function of wet size. The same optics are used
for MAM in GEOS-Chem.

The box is treated as a layer about 1 km thick (100 hPa), so the AOD values
correspond to a 1 km layer with the chosen aerosol load.

The optics use lookup tables read from NetCDF files (refractive indices per
species, Chebyshev fits per mode). These files must be present for the model
to run.

### Time evolution of temperature and humidity

Temperature and RH change **linearly** from their initial value at the first
step to their final value at the last step (namelist `mtmin`/`mrhmin` and
`mtmax`/`mrhmax`; the final value may be lower, e.g. for cooling). Set initial
= final for constant conditions; use different values to simulate e.g. a
humidification or cooling event.

---

## 2. Using the dashboard

### Starting it

```bash
cd GCMAMBOX/dashboard
python mam_dashboard.py
```

Then open <http://localhost:8050> in a browser.

To use another port (e.g. if 8050 is already taken, or to run two dashboards
side by side), pass `--port` (or `-p`):

```bash
python mam_dashboard.py --port 8051
```

and open <http://localhost:8051>.

### Setting up a case

The left panel, from top to bottom:

**Simulation Controls** — the **Run Simulation** button and a status line.
You can set the parameters manually, or start from a predefined case.

**Predefined Case** — loads an illustrative test case: it prefills every
parameter section below with the case's values. You can then change any value
and run. To reload a case after changing its values, select another case and
then the case again.

The cases go from one process at a time to more chemistry (nucleation,
SOA, black-carbon aging, ammonium nitrate, sea salt, dust, nitrate mode
shift). Each one is described briefly in the dashboard when selected.

**Time Parameters**
- *Time step* (s) and *Number of steps*. Total simulated time = step × number
  of steps (e.g. 300 s × 200 = 16.7 h). The default and the predefined cases
  use 300 s, the usual GEOS-Chem chemistry time step.
- *MOSAIC substeps* (1–10; 1 = standard MAM): splits each time step
  into that many substeps for the aerosol microphysics (gas–aerosol exchange
  with MOSAIC, renaming, nucleation, coagulation, aging), each of length
  time step / substeps. Gas production and the meteorology are still updated
  once per time step. In MAM, the NH3 produced during a step is given to
  MOSAIC all at once at the start of the step (H2SO4 is spread over the
  step); with substeps it is given in smaller portions. Use more substeps
  when gas or aerosol curves jump from one step to the next; they cost
  less than a shorter time step.
- At most 2000 steps per run, to keep runs short on a shared server. Long
  runs are animated with at most 150 frames (every k-th step), while the gas
  curves keep every step.

**Aerosol Processes**
- Tick the processes to include (see the table above). Switching processes
  off one at a time is the easiest way to see what each one does.

**Meteorological Parameters**
- *Temperature initial / final* in K, *Pressure* in Pa.
- *Relative humidity initial / final* as a **fraction between 0 and 1**
  (0.6 = 60 %), not in percent.

**Gas Phase**
- *Initial concentrations* (ppb) of H2SO4, HNO3, NH3, SOAG and HCl.
- *Production rates* (ppb/h, collapsible): constant sources of the same
  gases, applied at every step (0 = none). Without production the box is
  closed: the gases only partition between gas and particles. H2SO4 condenses
  within minutes, so a steady H2SO4 production rate (typically 0.01–0.1 ppb/h)
  is what drives nucleation and growth.

**Aerosol Initial Conditions**
- *Mode total mass* (µg/m³): initial dry aerosol mass in each of the 4 modes.
- *Mass fractions by mode*: composition of each mode. The fractions of a mode
  should add up to 1.
- The initial **number** of particles is not an input: it is derived from the
  mass, the composition (density) and the mode's default diameter.

**Accommodation Coefficients** (collapsible)
- Sticking probability (0–1) of H2SO4, HNO3, HCl and NH3 on particles. They
  control how fast these gases condense.

### Running

Click **Run Simulation**. The status line under the button reports success
or the error message from the model. A typical run takes well under a second;
runs are stopped after 30 s.

### Reading the results

The figure has four panels, animated over the time steps. The legend of each
panel gives the current value of each item (gas concentration, species mass,
mode number):

1. **Gas concentrations** vs time in hours (ppb, log scale).
2. **Aerosol mass size distribution** dm/dlogD (µg/m³) versus diameter (µm),
   stacked by species.
3. **Aerosol number size distribution** dN/dlogD (#/cm³) versus diameter,
   one curve per mode.
4. **Mode properties**: AOD at 550 nm per mode (wet aerosol, 1 km layer),
   stacked by species, with each mode's hygroscopicity, single scattering
   albedo, asymmetry parameter, longwave absorption and pH. pH is shown as
   N/A at the first step (initial state, before MOSAIC runs) and for modes
   without aerosol water.

Use **Play / Pause** to animate, or drag the slider to a given time step.
Hover over the curves to read values. To zoom, drag a box or drag the end of
an axis (e.g. extend the gas axis below 0.1 ppb to see H2SO4); double-click
to reset.

### Tips

- Start from a predefined case, change one parameter at a time, and compare.
- High RH increases water uptake, wet size and AOD, and favours nitrate and
  ammonium partitioning into the particles.
- An H2SO4 production rate drives new particle formation and growth of the
  Aitken mode, followed by renaming into the accumulation mode; NH3 and HNO3
  production rates feed ammonium nitrate formation.
- If a run fails, check that the mass fractions are sensible and that RH is
  given as a fraction.

### Adding a predefined case

Each case is a namelist file in `dashboard/assets/presets/`, which can also be
run directly with `gcmambox` (copy it to `namelist` in the run directory). To
add a case, put a new `.nml` file there starting with two comment lines,
`! title: ...` and `! description: ...`; after restarting the dashboard it
appears in the menu, ordered by file name.

---

## 3. Files

| Item | Location |
|---|---|
| Model source (Fortran) | `src/` (MAM library in `src/MAM/`) |
| Executable | `modbin/gcmambox` (build: `cd build && make`) |
| Optics tables | `opt_data/` (copies in `dashboard/assets/`) |
| Dashboard | `dashboard/mam_dashboard.py` |
| Predefined cases | `dashboard/assets/presets/*.nml` |
| Model output | `mam_output.nc` (NetCDF, written by each run) |

---

## 4. References

- Ghan, S. J. and Zaveri, R. A.: Parameterization of optical properties for
  hydrated internally mixed aerosol, *J. Geophys. Res.*, 112, D10201,
  <https://doi.org/10.1029/2006JD007927>, 2007.
- Liu, X., Easter, R. C., Ghan, S. J., Zaveri, R., et al.: Toward a minimal
  representation of aerosols in climate models: description and evaluation in
  the Community Atmosphere Model CAM5, *Geosci. Model Dev.*, 5, 709–739,
  <https://doi.org/10.5194/gmd-5-709-2012>, 2012.
- Wu, M., Wang, H., Easter, R. C., Lu, Z., Liu, X., Singh, B., et al.:
  Development and evaluation of E3SM-MOSAIC: Spatial distributions and
  radiative effects of nitrate aerosol, *J. Adv. Model. Earth Syst.*, 14,
  e2022MS003157, <https://doi.org/10.1029/2022MS003157>, 2022.
- Zaveri, R. A., Easter, R. C., Fast, J. D., and Peters, L. K.: Model for
  Simulating Aerosol Interactions and Chemistry (MOSAIC), *J. Geophys. Res.*,
  113, D13204, <https://doi.org/10.1029/2007JD008782>, 2008.
