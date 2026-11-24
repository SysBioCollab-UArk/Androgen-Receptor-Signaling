"""
ADT (Androgen Deprivation Therapy) Simulation Experiment
----------------------------------------------------------
Simulates the effect of an AR inhibitor (ARi) that competes with DHT
for binding at the AR ligand site. ARi represents enzalutamide-class
competitive AR antagonists used in clinical ADT.

ARi is added to the model as a new monomer in this script without
modifying AR_model.py. It binds competitively at the AR lig site,
preventing DHT binding and AR activation.

Conditions tested (ARi concentration at t=0, alongside 10 nM DHT):
  - Control: 0 nM ARi (no treatment)
  - 2.5 nM ARi
  - 5 nM ARi
  - 10 nM ARi  (1:1 equimolar with DHT)
  - 50 nM ARi
  - 100 nM ARi

Observables: Her2_p_tot, cPAcP_tot, ERK_pp_tot, AR_p_tot, PSA_tot, CycD_tot
Time window:  10-day DHT stimulation
"""

import os
import sys
import csv as _csv
import numpy as np
import glob

HOME = os.path.expanduser('~')
sys.path.insert(0, os.path.join(HOME, 'TIBD_microenvironment'))

PYDREAM_DIR = os.path.join(HOME, 'Androgen-Receptor-Signaling', 'Tasseff_2010')
sys.path.insert(0, PYDREAM_DIR)

# ── 1. Import model and add observables BEFORE importing run_ar_signaling_pydream
from AR_model import model
from pysb import Observable, Monomer, Parameter, Rule, Initial, ANY
from pysb.simulator import ScipyOdeSimulator

AR_mon  = model.monomers['AR']
ERK_mon = model.monomers['ERK']
CycD_mon = model.monomers['CycD']
Observable('AR_p_tot',   AR_mon(state='p'))
Observable('ERK_pp_tot', ERK_mon(state='pp'))
Observable('CycD_tot',   CycD_mon())

# ── 2. Add ARi monomer and competitive binding rules to the model
#       ARi binds at the AR lig site (same site as DHT), preventing DHT binding.
#       Binding kinetics initialized to DHT values as a baseline assumption.
#       These parameters are NOT calibrated — they represent fixed drug properties.
#       Components are defined at module scope (same pattern as Observable additions)
#       so PySB's SelfExporter automatically attaches them to the active model.
Monomer('ARi', ['b'])
Parameter('kf_AR_binds_ARi', 2.486)     # same as kf_AR_binds_DHT (baseline assumption)
Parameter('kr_AR_binds_ARi', 5.42e-5)   # same as kr_AR_binds_DHT (baseline assumption)
Parameter('ARi_0', 0)                    # start at zero; injected via protocol
ARi = model.monomers['ARi']
Rule('AR_binds_ARi',
     AR_mon(lig=None) + ARi(b=None) | AR_mon(lig=1) % ARi(b=1),
     model.parameters['kf_AR_binds_ARi'],
     model.parameters['kr_AR_binds_ARi'])
Initial(ARi(b=None), model.parameters['ARi_0'])

# ── 3. Import calibration infrastructure
from run_ar_signaling_pydream import (expt_data_file, sim_protocols,
                                      custom_priors, no_sample, param_expts_map,
                                      observables as obs_list)
from param_calibration import ParameterCalibration
from SIM_PROTOCOLS.sim_protocols import SequentialInjections, ScaleBkProtocol

# Exclude ARi parameters from Bayesian sampling — they are fixed drug properties,
# not model parameters to be inferred from the calibration data.
no_sample = list(no_sample) + ['ARi_0', 'kf_AR_binds_ARi', 'kr_AR_binds_ARi']

calibrator = ParameterCalibration(model,
                                  expt_data_file,
                                  sim_protocols,
                                  priors=custom_priors,
                                  no_sample=no_sample,
                                  param_expts_map=param_expts_map)

# ── 4. Load and filter posterior samples
sample_files = sorted(glob.glob(os.path.join(PYDREAM_DIR,
                                             'dreamzs_5chain_sampled_params_chain_*.npy')))
print(f'Found {len(sample_files)} sampled_params files')

param_samples, log_ps = calibrator.plot_param_dist(
    sample_files,
    cutoff=2,
    show_plot=False,
    save_plot=False
)
print(f'Filtered posterior samples shape: {param_samples.shape}')

MAX_SAMPLES = 2000
if param_samples.shape[0] > MAX_SAMPLES:
    rng = np.random.default_rng(42)
    idx = rng.choice(param_samples.shape[0], size=MAX_SAMPLES, replace=False)
    param_samples = param_samples[idx]
    log_ps = log_ps[idx]
    print(f'Subsampled to {MAX_SAMPLES} posterior samples for simulation')

base_idxs = list(calibrator.samples_idxs[0])

# ── 5. Build simulation protocols for each ARi dose
#       Each condition uses the same calibrated parameters but a different
#       ARi injection amount at t=0 alongside the standard 10 nM DHT.
#       expt_data is pulled from the existing calibration protocol so
#       ScaleBkProtocol can compute Bk scaling on the same scale as calibration.
t_equil = 10 * 24 * 3600
ARi_doses_nM = [0, 2.5, 5, 10, 50, 100]
condition_labels = [f'{d} nM ARi' if d > 0 else 'control (0 nM ARi)'
                    for d in ARi_doses_nM]

expt_data = sim_protocols[0].expt_data

solver = ScipyOdeSimulator(model)

adt_protocols = []
for dose in ARi_doses_nM:
    inj = SequentialInjections(
        solver,
        t_equil=t_equil,
        time_perturb_value={0: [('DHT(b=None)', 10), ('ARi(b=None)', dose)]}
    )
    adt_protocols.append(ScaleBkProtocol(inj, obs_list, expt_data=expt_data))

# ── 6. Run simulations
tspan_days = np.linspace(0, 10 * 24 * 3600, 101)
observables = ['Her2_p_tot', 'cPAcP_tot', 'ERK_pp_tot', 'AR_p_tot', 'PSA_tot', 'CycD_tot']
SAVE_DIR = os.path.dirname(__file__)
n_cond = len(adt_protocols)

def _save_with_mean(protocols, tspan, samples, parameter_idxs,
                    samples_idxs_list, observables, save_path):
    base_pv = np.array([p.value for p in model.parameters])
    if not isinstance(protocols, list):
        protocols = [protocols] * len(samples_idxs_list)
    rows = []
    for sim_id, (protocol, idxs) in enumerate(zip(protocols, samples_idxs_list)):
        print(f'\nMean pass: condition {sim_id + 1}/{len(samples_idxs_list)}')
        yvals_dict = {obs: [] for obs in observables}
        n_removed = 0
        for i, sample in enumerate(samples):
            pv = base_pv.copy()
            pv[parameter_idxs] = 10 ** sample[idxs]
            output = protocol.run(tspan, pv)
            if np.any(np.isnan(output[observables[0]])):
                n_removed += 1
                continue
            for obs in observables:
                yvals_dict[obs].append(output[obs])
            if (i + 1) % 200 == 0:
                print(i + 1, end=' ', flush=True)
        print(f' done ({n_removed} NaN outputs skipped)')
        for obs in observables:
            yvals = np.array(yvals_dict[obs])
            ymin  = np.percentile(yvals, 5,  axis=0)
            ymax  = np.percentile(yvals, 95, axis=0)
            ymean = np.mean(yvals, axis=0)
            for t, lo, hi, mu in zip(tspan / 86400, ymin, ymax, ymean):
                rows.append([obs, t, lo, hi, mu, sim_id])
    with open(f'{save_path}.csv', 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(['observable', 'time', 'yval_min', 'yval_max', 'yval_mean', 'sim_id'])
        w.writerows(rows)
    print(f'Saved: {save_path}.csv')

print('\n=== Running ADT simulations ===')
calibrator.plot_timecourses(
    model,
    tspans         = [tspan_days] * n_cond,
    sim_protocols  = adt_protocols,
    param_samples  = param_samples,
    parameter_idxs = calibrator.parameter_idxs,
    observables    = [observables] * n_cond,
    samples_idxs   = [base_idxs] * n_cond,
    exp_data       = None,
    show_plot      = False,
    save_plot      = SAVE_DIR,
    separate_plots = False,
    xlabel         = 'time (s)',
    ylabels        = [['amount (relative)'] * len(observables)] * n_cond,
    leg_labels     = [[label] * len(observables) for label in condition_labels],
    save_sim_data  = False
)
print('ADT simulations DONE')
_save_with_mean(
    adt_protocols, tspan_days, param_samples, calibrator.parameter_idxs,
    [base_idxs] * n_cond, observables,
    os.path.join(SAVE_DIR, 'SIM_DATA_ADT')
)
print('\nAll simulations complete. Plots saved to:', SAVE_DIR)
