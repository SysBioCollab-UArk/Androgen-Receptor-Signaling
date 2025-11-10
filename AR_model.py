from pysb import *
from util import set_model, create_transcription_rules, create_translation_rules

Model()
set_model(model)  # put the model into the global namespace

# monomers

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
Monomer('ERK',    ['mek','pase3','sos','ets','ap1','state'],{'state':['u','p','pp']})
Monomer('Pase2',  ['mek'])
Monomer('Pase3',  ['erk'])

# 4) ERK nuclear targets
Monomer('ETS',    ['erk','gene','state'],{'state':['u','p']})
Monomer('AP1',    ['erk','gene','state'],{'state':['u','p']})

# 5) Internalization forms used explicitly in reactions
# Monomer('EGFRi',  ['l','d','grb2','sos','shc','erk'])
# Monomer('EGFi')

# 6) PI3K–AKT–mTOR arm
Monomer('PI3K',   ['egfr_her2','ptdins2','sos','state'],{'state':['i','act']})
Monomer('PtdIns2',['pi3k'])  # PIP2
Monomer('PtdIns3',['pten','akt','pdk1'])  # PIP3
Monomer('PTEN',   ['ptdins3'])
Monomer('Akt',    ['ptdins3','pdk1','tor','pase7','state'],{'state':['i','m','act']})
Monomer('Pdk1',   ['ptdins3','akt','state'],{'state':['i','m']})
Monomer('TOR',    ['akt','e4ebp1','state'],{'state':['i','act']})
Monomer('_4EBP1', ['eif4e','tor','state'],{'state':['u','p']})
Monomer('Pase7',  ['akt'])

# 7) Transport/secreted forms used explicitly
# Monomer('T_e')
# Monomer('sPAcP_e')

# 8) AR axis and modifiers
Monomer('AR',     ['hsp','t','DHT','erk','ar','pase5','state'],{'state':['u','p']})
Monomer('T',      ['rase5a','ar','loc'], {'loc': ['intra', 'extra']})
Monomer('Rase5a', ['t'])
Monomer('HSP',    ['ar'])
Monomer('DHT',    ['ar'])
Monomer('Pase5',  ['ar','ets','ap1'])
Monomer('Pase6',  ['ap1'])

# 9) Transcription (promoters/genes) & RNA Pol II
Monomer('RNAp',   ['gene'])
Monomer('g_cPAcP',['rnap','tf'])
Monomer('g_sPAcP',['rnap','tf'])
Monomer('g_CycD', ['rnap','tf'])
Monomer('g_PSA',  ['rnap','tf'])

# 10) mRNAs and translation machinery
Monomer('eIF4E',     ['mrna_4ebp1'])
Monomer('_40S',      ['mrna', '_60s'])
Monomer('_60S',      ['_40s'])
Monomer('mRNA_cPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_sPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_CycD', ['eif4e','_40s','elong'], {'elong': ['i','a']})
Monomer('mRNA_PSA',  ['eif4e','_40s','elong'], {'elong': ['i','a']})

# 11) Protein products (sinks in your list)
Monomer('CycD')
Monomer('PSA')


# initial concentrations
Parameter('Lig_0', 100)
Parameter('Rec_0', 100)

Initial(EGF(r=None, loc='mem'), Lig_0)
Initial(EGFR(l=None, d=None, grb2_shc=None, state='u', loc='mem'), Rec_0)


# rules

# CycD transcription
tfs = [ETS(erk=None,gene=None,state='p'), AP1(erk=None,gene=None,state='p')]
kf_kr_kcat = [
    [4.952e-5, 0.1104, 1.099e-2],  # basal
    [0.2931, 9.037e-3, 1.156e-2],  # ETS-p
    [0.4288, 2.945e-2, 3.292e-2]]  # AP1-p
k_tf_on_off = [[0.138, 1.26], [0.3726, 2.171]]  # ETS-p, AP1-p
k_deg = 0.8094
create_transcription_rules(CycD, kf_kr_kcat, k_deg, tfs=tfs, k_tf_on_off=k_tf_on_off)

# CycD translation
kf_kr = [
    [2.137e-2, 8.14e-3],   # eIF4E
    [8.773e-2, 5.877e-3],  # 40S
    [0.7418, 1.306e-3]     # 60S
]
k_release = 1.194
k_elongate = 0.9399
k_terminate = 1.781
k_deg = 6.123e-3
create_translation_rules(CycD, kf_kr, k_release, k_elongate, k_terminate, k_deg)


# Add other transcription and translation rules
# ...


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
