from AR_model import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *
import os

this_dir = os.path.dirname(__file__)
expt_data_file = os.path.join(this_dir, 'DATA', 'Tasseff_2010.csv')
expt_data = pd.read_csv(expt_data_file)

solver = ScipyOdeSimulator(model)

DHT_stimulation_10_nM = SequentialInjections(solver, t_equil=10 * 24 * 3600,  # 10 days
                                                 time_perturb_value={0: ('DHT(b=None)', 10)})
observables = ['Her2_p_tot', 'cPAcP_tot', 'PSA_tot']
protocol_A = ScaleBkProtocol(DHT_stimulation_10_nM, observables, expt_data=expt_data)

custom_priors = None

no_sample = [
    # initial conditions
    'EGF_loc_extra_0', 'EGF_loc_intra_0', 'EGFR_state_u_loc_intra_0', 'EGFR_state_p_loc_extra_0',
    'EGFR_state_p_loc_intra_0', 'Her2_state_p_0', 'Ras_state_GTP_0', 'Shc_state_p_0', 'cPAcP_0', 'sPAcP_loc_intra_0',
    'sPAcP_loc_extra_0', 'Raf_state_p_0', 'MEK_state_p_0', 'MEK_state_pp_0', 'ERK_state_p_0', 'ERK_state_pp_0',
    'ETS_state_p_0', 'AP1_state_p_0', 'PI3K_state_act_0', 'Akt_state_m_0', 'Akt_state_act_0', 'Pdk1_state_m_0',
    'TOR_state_act_0', '_4EBP1_state_p_0', 'AR_state_p_0', 'T_loc_intra_0', 'T_loc_extra_0', 'DHT_0',
    'mRNA_cPAcP_elong_i_0', 'mRNA_cPAcP_elong_a_0', 'mRNA_sPAcP_elong_i_0', 'mRNA_sPAcP_elong_a_0',
    'mRNA_CycD_elong_i_0', 'mRNA_CycD_elong_a_0', 'mRNA_PSA_elong_i_0', 'mRNA_PSA_elong_a_0', 'CycD_0', 'PSA_0',
    # rate constants
    'k_EGF_internalize']

sim_protocols = [protocol_A]
param_expts_map=None


if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      expt_data_file,
                                      sim_protocols,
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)

    calibrator.run(niterations=50000, nchains=5, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True}, restart=False)
