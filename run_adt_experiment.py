"""
ADT (Androgen Deprivation Therapy) Simulation Experiment
----------------------------------------------------------
Simulates the effect of an AR inhibitor (ARi) that competes with DHT for binding at the AR ligand site. ARi represents
enzalutamide-class competitive AR antagonists used in clinical ADT.

ARi is added to the model as a new monomer in this script without modifying AR_model.py. It binds competitively at the
AR lig site, preventing DHT binding and AR activation.

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

THIS_DIR = os.path.dirname(__file__)
PYDREAM_DIR = os.path.join(THIS_DIR, 'PyDREAM_Output')
sys.path.insert(0, PYDREAM_DIR)

from run_ar_signaling_pydream import *
from pysb import Observable, Monomer, Parameter, Rule, Initial
from pysb.util import alias_model_components
import csv
import numpy as np

model.reset_equations()  # clear out model.species, model.reactions, etc.

# ── 1. Add observables to model imported from `run_ar_signaling_pydream`
alias_model_components()

Observable('AR_p_tot',   AR(state='p'))
Observable('ERK_pp_tot', ERK(state='pp'))
Observable('CycD_tot',   CycD())

# ── 2. Add ARi monomer and competitive binding rules to the model.
#       ARi binds at the AR lig site (same site as DHT), preventing DHT binding. Binding kinetics initialized to DHT
#       values as a baseline assumption. These parameters are NOT calibrated — they represent fixed drug properties.
#       Components are defined at module scope (same pattern as Observable additions), so PySB's SelfExporter
#       automatically attaches them to the active model.
Monomer('ARi', ['b'])
Parameter('kf_AR_binds_ARi', 2.486)     # same as kf_AR_binds_DHT (baseline assumption)
Parameter('kr_AR_binds_ARi', 5.42e-5)   # same as kr_AR_binds_DHT (baseline assumption)
Parameter('ARi_0', 0)                   # start at zero; injected via protocol

alias_model_components()

Rule('AR_binds_ARi', AR(lig=None) + ARi(b=None) | AR(lig=1) % ARi(b=1), kf_AR_binds_ARi, kr_AR_binds_ARi)
Initial(ARi(b=None), ARi_0)

# ── 3. Build simulation protocols for each ARi dose.
#       Each condition uses the same calibrated parameters but a different ARi injection amount at t=0 alongside the
#       standard 10 nM DHT. `expt_data` is imported from `run_ar_signaling_pydream`.
ARi_doses_nM = [0, 2.5, 5, 10, 50, 100]
condition_labels = [f'{d} nM ARi' if d > 0 else 'control (0 nM ARi)' for d in ARi_doses_nM]

t_equil = 10 * 24 * 3600  # 10 days
solver = ScipyOdeSimulator(model)

adt_protocols = []
for dose in ARi_doses_nM:
    inj = SequentialInjections(solver, t_equil=t_equil,
                               time_perturb_value={0: [('DHT(b=None)', 10), ('ARi(b=None)', dose)]})
    adt_protocols.append(ScaleBkProtocol(inj, observables, expt_data=expt_data))

# ── 4. Create the calibrator
no_sample += ['kf_AR_binds_ARi', 'kr_AR_binds_ARi', 'ARi_0']
calibrator = ParameterCalibration(model,
                                  expt_data_file,
                                  sim_protocols=adt_protocols,
                                  priors=custom_priors,
                                  no_sample=no_sample,
                                  param_expts_map=param_expts_map)

# ── 5. Load and filter posterior samples
kwargs = {'show_plots': True,
          'plot_ll_args': {'cutoff': 2, 'fontsizes': {'labels': 14, 'ticks': 14, 'legend': 14}, 'leg_ncols': 2},
          'plot_pd_args': {'sharex': 'all', 'fontsize': 28, 'labelsize': 24, 'bw_adjust': 3, 'plot_priors': True},
          'plot_tc_args': {'separate_plots': False, 'save_sim_data': True},
          'which_plots': 2}

max_iter = 50000  # Change this to `None` to use all the PyDREAM samples
_, (param_samples, _, _, _), _ = plot_pydream_output(PYDREAM_DIR, calibrator, max_iter=max_iter, **kwargs)

# ── 6. Run simulations
tspan_days = np.linspace(0, 10 * 24 * 3600, 101)
observables += ['ERK_pp_tot', 'AR_p_tot', 'CycD_tot']  # First three are ['Her2_p_tot', 'cPAcP_tot', 'PSA_tot']
n_cond = len(adt_protocols)

ParameterCalibration.plot_timecourses(
    model,
    tspans         = [tspan_days] * n_cond,
    sim_protocols  = adt_protocols,
    param_samples  = param_samples,
    parameter_idxs = calibrator.parameter_idxs,
    observables    = [observables] * n_cond,
    samples_idxs   = [calibrator.samples_idxs[0]] * n_cond,
    exp_data       = calibrator.raw_data,
    show_plot      = True,
    save_plot      = PYDREAM_DIR,
    separate_plots = False,
    xlabel         = 'time (s)',
    ylabels        = [['amount'] * len(observables)] * n_cond,
    leg_labels     = [[label] * len(observables) for label in condition_labels],
    save_sim_data  = os.path.join(PYDREAM_DIR, 'SIM_DATA'),
    num_tc         = 10
)
