from pysb import *
from util import set_model, create_transcription_rules, create_translation_rules

Model()
set_model(model)  # put the model into the global namespace

# === MONOMERS ===

# 1) Ligands & Receptors
Monomer('EGF',  ['r','loc'], {'loc': ['mem','cyt']})
Monomer('EGFR', ['l','d','grb2_shc','state','loc'],{'state':['u','p'], 'loc': ['mem','cyt']})
Monomer('Her2', ['d','grb2_shc','cpacp','state'],{'state':['u','p']})
Monomer('Grb2', ['sos','shc','egfr_her2'])
Monomer('Sos',  ['grb2','ras_erk','pi3k'])
Monomer('Ras',  ['sos','gap','raf','state'],{'state':['GDP','GTP']})
Monomer('Shc',  ['egfr_her2','grb2','state'],{'state':['u','p']})

# 2) Phosphatases acting on HER2/AR modules
Monomer('cPAcP',  ['her2']) #dimer
Monomer('sPAcP',  ['her2','loc'], {'loc': ['intra','extra']})

# 3) RAS–RAF–MEK–ERK cascade
Monomer('GAP',    ['ras'])
Monomer('Raf',    ['ras','mek','pase1','state'],{'state':['u','p']})
Monomer('Pase1',  ['raf'])
Monomer('MEK',    ['raf','erk','pase2','state'],{'state':['u','p','pp']})
Monomer('ERK',    ['ar', 'mek','pase3','sos','ets','ap1','state'],{'state':['u','p','pp']})
Monomer('Pase2',  ['mek'])
Monomer('Pase3',  ['erk'])

# 4) ERK nuclear targets
Monomer('ETS',    ['erk','gene','state'],{'state':['u','p']})
Monomer('AP1',    ['erk','gene','state'],{'state':['u','p']})

# 5) PI3K–AKT–mTOR arm
Monomer('PI3K',   ['egfr_her2','ptdins2','sos','state'],{'state':['i','act']})
Monomer('PtdIns2',['pi3k'])  # PIP2
Monomer('PtdIns3',['pten','akt','pdk1'])  # PIP3
Monomer('PTEN',   ['ptdins3'])
Monomer('Akt',    ['ptdins3','pdk1','tor','pase7','state'],{'state':['i','m','act']})
Monomer('Pdk1',   ['ptdins3','akt','state'],{'state':['i','m']})
Monomer('TOR',    ['akt','e4ebp1','state'],{'state':['i','act']})
Monomer('_4EBP1', ['eif4e','tor','state'],{'state':['u','p']})
Monomer('Pase7',  ['akt'])

# 6) AR axis and modifiers
Monomer('AR',     ['lig','ar','erk','pase5','gene','state'], {'state':['u','p']})
Monomer('T',      ['b','loc'], {'loc': ['intra', 'extra']})  # testosterone
Monomer('Rase5a', ['t'])
Monomer('HSP',    ['ar'])
Monomer('DHT',    ['b'])  # dihydrotestosterone
Monomer('Pase5',  ['ar','ets','ap1'])
Monomer('Pase6',  ['ap1'])

# 7) Transcription (promoters/genes) & RNA Pol II
Monomer('RNAp',   ['gene'])
Monomer('g_cPAcP',['rnap','tf1','tf2'])
Monomer('g_sPAcP',['rnap','tf'])
Monomer('g_CycD', ['rnap','tf'])
Monomer('g_PSA',  ['rnap','tf'])

# 8) mRNAs and translation machinery
Monomer('eIF4E',     ['mrna_4ebp1'])
Monomer('_40S',      ['mrna', '_60s'])
Monomer('_60S',      ['_40s'])
Monomer('mRNA_cPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_sPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_CycD', ['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_PSA',  ['eif4e','_40s','elong'], {'elong': ['i','a']})

# 9) Protein products (sinks in your list)
Monomer('CycD')
Monomer('PSA')

# === INITIAL CONCENTRATIONS ===

Parameter('Lig_0', 100)
Parameter('Rec_0', 100)

Initial(EGF(r=None, loc='mem'), Lig_0)
Initial(EGFR(l=None, d=None, grb2_shc=None, state='u', loc='mem'), Rec_0)

# === RULES ===

# CycD transcription
kf_kr_kcat = [[4.952e-5, 0.1104, 1.099e-2]]  # basal
k_mRNA_deg = 0.8094
tf_species = [ETS(erk=None,gene=None,state='p'), AP1(erk=None,gene=None,state='p')]
tf_names = ['ETS_p', 'AP1_p']
k_tf_on_off = [[0.138, 1.26], [0.3726, 2.171]]  # ETS-p, AP1-p
kf_kr_kcat.extend(
    [[0.2931, 9.037e-3, 1.156e-2],  # ETS-p
    [0.4288, 2.945e-2, 3.292e-2]]   # AP1-p
)
create_transcription_rules(CycD, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species, tf_names)), k_tf_on_off=k_tf_on_off)

# CycD translation
kf_kr = [
    [2.137e-2, 8.14e-3],   # eIF4E
    [8.773e-2, 5.877e-3],  # 40S
    [0.7418, 1.306e-3]     # 60S
]
k_release = 1.194
k_elongate = 0.9399
k_terminate = 1.781
k_prot_deg = 6.123e-3
create_translation_rules(CycD, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)

# PSA transcription
# ... TODO

# PSA translation
# ... TODO

# cPAcP transcription
kf_kr_kcat = [[5.788e-2, 7.556e-2, 1.311e-2]]  # basal
k_mRNA_deg = 1.3
tf_species = []
tf_names = []
k_tf_on_off = []
# AR-p-2 + g-cPAcP <-> AR-p-2-g-cPAcP
tf_species.append(
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_2')
k_tf_on_off.append([1.618e-3, 9.872e-5])
# AR-p-DHT-2 + g-cPAcP <-> AR-p-DHT-2-g-cPAcP
tf_species.append(
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    DHT(b=3) % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_DHT_2')
k_tf_on_off.append([8.423e-4, 1.366e-7])
# AR-p-DHT-AR-p + g-cPAcP <-> AR-p-DHT-AR-p-g-cPAcP
tf_species.append(
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_DHT_AR_p')
k_tf_on_off.append([4.544e-3, 1.314e-7])
# AR-p-T-2 + g-cPAcP <-> AR-p-T-2-g-cPAcP
tf_species.append(
    T(b=2, rase5a=None, loc='intra') % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    T(b=3, rase5a=None, loc='intra') % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_T_2')
k_tf_on_off.append([6.374e-2, 2.954e-4])
# AR-p-T-AR-p + g-cPAcP <-> AR-p-T-AR-p-g-cPAcP
tf_species.append(
    T(b=2, rase5a=None, loc='intra') % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_T_AR_p')
k_tf_on_off.append([6.191e-2, 1.099e-4])
# AR-p-DHT-AR-p-T + g-cPAcP <-> AR-p-DHT-AR-p-T-g-cPAcP
tf_species.append(
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    T(b=3, rase5a=None, loc='intra') % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p'))
tf_names.append('AR_p_DHT_AR_p_T')
k_tf_on_off.append([0.1543,   2.343e-4])
create_transcription_rules(cPAcP, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species, tf_names)), k_tf_on_off=k_tf_on_off)

# cPAcP translation
kf_kr = [
    [1.351e-2, 9.892e-3],  # eIF4E
    [0.233, 9.908e-3],     # 40S
    [0.9053, 1.499e-3]     # 60S
]
k_release = 1.185
k_elongate = 1.31
k_terminate = 1.612
k_prot_deg = 1.208e-2
create_translation_rules(cPAcP, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)

# sPAcP transcription
# ... TODO

# sPAcP translation
# ... TODO

# AR+ERK-pp	↔	AR-ERK-pp
# AR-ERK-pp	↔	AR-p+ERK-pp
Parameter('kf_ERKpp_phos_AR', 1.873e-3)
Parameter('kr_AERKpp_phos_AR', 0.388)
Parameter('kcat_ERKpp_phos_AR', 2.57E-2)
Rule('AR_binds_ERKpp',
     AR(lig=None, ar=None, erk=None, pase5=None, state='u') +
     ERK(ar=None, mek=None, pase3=None, sos=None, ets=None, ap1=None, state='pp') |
     AR(lig=None, ar=None, erk=1, pase5=None, state='u') %
     ERK(ar=1, mek=None, pase3=None, sos=None, ets=None, ap1=None, state='pp'),
     kf_ERKpp_phos_AR, kr_AERKpp_phos_AR)
Rule('AR_ERKpp_phosphorylates_AR',
     AR(lig=None, ar=None, erk=1, pase5=None, state='u') %
     ERK(ar=1, mek=None, pase3=None, sos=None, ets=None, ap1=None, state='pp') >>
     AR(lig=None, ar=None, erk=None, pase5=None, state='p') +
     ERK(ar=None, mek=None, pase3=None, sos=None, ets=None, ap1=None, state='pp'),
     kcat_ERKpp_phos_AR)

# T+Rase5a	↔	T-Rase5a
# T-Rase5a	→	DHT+Rase5a
Parameter('kf_T_Rase5a', 2.171e-2)
Parameter('kr_T_Rase5a', 0.3308)
Parameter('kcat_T_Rase5a', 1.374)
Rule('T_binds_Rase5a',
     T(b=None, loc='intra') + Rase5a(t=None) | T(b=1, loc='intra') % Rase5a(t=1), kf_T_Rase5a, kr_T_Rase5a)
Rule('T_converts_to_DHT', T(b=1, loc='intra') % Rase5a(t=1) >> DHT(b=None) + Rase5a(t=None), kcat_T_Rase5a)

# AR+HSP	↔	AR-HSP
Parameter('kf_AR_binds_HSP', 9.162e-3)
Parameter('kr_AR_binds_HSP', 1.01e-3)
Rule('AR_binds_HSP',
     AR(lig=None, ar=None, erk=None, pase5=None, state='u') + HSP(ar=None) |
     AR(lig=1, ar=None, erk=None, pase5=None, state='u') % HSP(ar=1),
     kf_AR_binds_HSP, kr_AR_binds_HSP)

# AR+T	↔	AR-T
# AR-T	→	AR-p-T
Parameter('kf_AR_binds_T', 5.289e-3)
Parameter('kr_AR_binds_T', 9.361e-4)
Parameter('kcat_AR_binds_T', 5.458E-2)
Rule('AR_binds_T',
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='u') + T(b=None) |
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='u') % T(b=1),
     kf_AR_binds_T, kr_AR_binds_T)
Rule('AR_T_phosphorylates',
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='u') % T(b=1) >>
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=1), kcat_AR_binds_T)

# AR-p+AR-p-T	↔	AR-p-AR-p-T
Parameter('kon_AR_p_AR_p_T', 0.532)
Parameter('koff_AR_p_AR_p_T', 6.128e-4)
Rule('AR_p_binds_AR_p_T',
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='p') +
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=1) |
     AR(lig=None, ar=2, erk=None, pase5=None, gene=None, state='p') +
     AR(lig=1, ar=2, erk=None, pase5=None, gene=None, state='p') % T(b=1),
     kon_AR_p_AR_p_T, koff_AR_p_AR_p_T)

# 2*AR-p-T	↔	AR-p-T-2
Parameter('kon_AR_p_T_dimerize', 0.5835)
Parameter('koff_AR_p_T_dimerize', 5.149e-4)
Rule('AR_p_T_dimerizes',
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=1) +
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=1) |
     AR(lig=1, ar=2, erk=None, pase5=None, gene=None, state='p') % T(b=1) +
     AR(lig=1, ar=2, erk=None, pase5=None, gene=None, state='p') % T(b=1),
     kon_AR_p_T_dimerize, koff_AR_p_T_dimerize)

# 2*AR-p	↔	AR-p-2
Parameter('kon_AR_p_dimerize', 0.2848)
Parameter('koff_AR_p_dimerize', 0.1281)
Rule('AR_p_dimerizes',
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='p') +
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='p') |
     AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p') +
     AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'),
     kon_AR_p_dimerize, koff_AR_p_dimerize)

# AR+DHT	↔	AR-DHT
# AR-DHT	→	AR-p-DHT

# AR-p-DHT+AR-p-T	↔	AR-p-DHT-AR-p-T

# AR-p-DHT+AR-p	↔	AR-p-DHT-AR-p

# 2*AR-p-DHT	↔	AR-p-DHT-2

# AR-p+Pase5	↔	AR-p-Pase5
# AR-p-Pase5	→	AR+Pase5

# AR-p-T+Pase5	↔	AR-p-T-Pase5
# AR-p-T-Pase5	→	AR-T+Pase5

# AR-p-DHT+Pase5	↔	AR-p-DHT-Pase5
# AR-p-DHT-Pase5	→	AR-DHT+Pase5

# T-e	↔	T
Parameter('kf_T_extra_to_intra', 1.449)
Parameter('kr_T_extra_to_intra', 1.555)
Rule('T_extra_to_intra', T(b=None, loc='extra') | T(b=None, loc='intra'), kf_T_extra_to_intra, kr_T_extra_intra)


Parameter('kf_LR', 1)
Parameter('kr_LR', 1)
Rule('lig_binds_rec', EGF(r=None) + EGFR(l=None, d=None) | EGF(r=1) % EGFR(l=1, d=None), kf_LR, kr_LR)

Parameter('kf_Her2_cPAcP', 1)
Parameter('kr_Her2_cPAcP', 1)
Rule('Her2p_binds_cPAcP',
     Her2(d=ANY, cpacp=None,state='p') + cPAcP(her2=None) | Her2(d=ANY, cpacp=1, state='p') % cPAcP(her2=1),
     kf_Her2_cPAcP, kr_Her2_cPAcP)

Parameter('kf_dimer', 1)
Parameter('kr_dimer', 1)
Rule('EGFR_dimerization', EGFR(l=ANY, d=None)+ EGFR(l=ANY, d=None) | EGFR(l=ANY, d=1) % EGFR(l=ANY, d=1),
     kf_dimer, kr_dimer)

Parameter('kf_phos', 1)
Rule('EGFR_autophosphorylation', EGFR(d=ANY, state='u') >> EGFR(d=ANY, state='p'), kf_phos)

Parameter('kf_EGFR_Shc', 1)
Parameter('kr_EGFR_Shc', 1)
Rule('EGFRp_binds_Shc',
     EGFR(state='p', grb2_shc=None) + Shc(egfr_her2=None) | EGFR(state='p', grb2_shc=1) % Shc(egfr_her2=1),
     kf_EGFR_Shc, kr_EGFR_Shc)

Parameter('kf_Shc_Grb2', 1)
Parameter('kr_Shc_Grb2', 1)
Rule('Shcp_binds_Grb2', Shc(state='p', grb2=None) + Grb2(shc=None) | Shc(state='p', grb2=1) % Grb2(shc=1),
     kf_Shc_Grb2, kr_Shc_Grb2)


# observables
Observable('Lig_free', EGF(r=None))
Observable('Lig_bound', EGF(r=ANY))
Observable('Rec_unphos', EGFR(state='u'))
Observable('Rec_phos', EGFR(state='p'))

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from pysb.simulator import ScipyOdeSimulator

    # simulation commands
    tspan = np.linspace(0, 0.5, 101)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    result = sim.run()

    # plot results
    plt.figure(constrained_layout=True)
    for obs in model.observables:
        plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.legend(loc='best')

    plt.show()
