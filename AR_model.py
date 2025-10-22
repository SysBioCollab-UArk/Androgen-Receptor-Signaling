from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

# monomers

# Ligands & Receptors
Monomer('EGF', ['r'])
Monomer('EGFR', ['l','d','grb2_shc','state'],{'state':['u','p']})
Monomer('Her2', ['d','grb2_shc','cpacp','state'],{'state':['u','p']})
Monomer('Grb2', ['sos','shc','egfr_her2'])
Monomer('Shc', ['egfr_her2','grb2','state'],{'state':['u','p']})
Monomer('Sos', ['grb2','ras_erk','pi3k'])
Monomer('Ras', ['sos','gap','raf','state'],{'state':['GDP','GTP']})


# 2) Phosphatases acting on HER2/AR modules
Monomer('cPAcP',  ['her2']) #dimer
Monomer('sPAcP',  ['her2'])

# 3) RAS–RAF–MEK–ERK cascade
Monomer('GAP',    ['ras'])
Monomer('Raf',    ['ras','mek','pase1','state'],{'state':['u','p']})
Monomer('Pase1',  ['raf'])
Monomer('MEK',    ['raf','erk','pase2','state'],{'state':['u','p','pp']})
Monomer('Pase2',  ['mek'])
Monomer('ERK',    ['mek','pase3','sos','ets','ap1','state'],{'state':['u','p','pp']})
Monomer('Pase3',  ['erk'])

# 4) ERK nuclear targets
Monomer('ETS',    ['erk','rnap','state'],{'state':['u','p']})
Monomer('AP1',    ['erk','rnap','state'],{'state':['u','p']})

# 5) Internalization forms used explicitly in reactions
Monomer('EGFRi',  ['l','d','grb2','sos','shc','erk'])
Monomer('EGFi')

# 6) PI3K–AKT–mTOR arm
Monomer('PI3K',   ['egfr_her2','ptdins2','sos','state'],{'state':['i','act']})
Monomer('PtdIns2',['pi3k'])
Monomer('PtdIns3',['pten','akt','pdk1'])
Monomer('PTEN',   ['ptdins3'])
Monomer('Akt',    ['ptdins3','pdk1','tor','pase7','state'],{'state':['i','m','act']})
Monomer('Pdk1',   ['ptdins3','akt','state'],{'state':['i','m']})
Monomer('TOR',    ['akt','e4ebp1','state'],{'state':['i','act']})
Monomer('E4EBP1', ['eif4e','tor','state'],{'state':['u','P']})
Monomer('eIF4E',  ['e4ebp1','mrna'])
Monomer('Pase7',  ['akt'])

# 7) Transport/secreted forms used explicitly
Monomer('T_e')
Monomer('T',      ['rase5a','ar'])
Monomer('sPAcP_e')

# 8) AR axis and modifiers
Monomer('Rase5a', ['T'])
Monomer('DHT',    ['ar'])
Monomer('AR',     ['hsp','T','DHT','erk','ar','pase5','state'],{'state':['u','p']})
Monomer('HSP',    ['ar'])
Monomer('Pase5',  ['ar','ets','ap1'])
Monomer('Pase6',  ['ap1'])

# 9) Transcription (promoters/genes) & RNA Pol II
Monomer('RNAp',   ['promoter'])
Monomer('g_CycD', ['rnap','ets_p','ap1_p'])
Monomer('g_PSA',  ['rnap','ar_dimer'])
Monomer('g_cPAcP',['rnap'])
Monomer('g_sPAcP',['rnap'])

# 10) mRNAs and translation machinery
Monomer('mRNA_CycD',['eif4e'])
Monomer('mRNA_PSA', ['eif4e'])
Monomer('mRNA_cPAcP',['eif4e'])
Monomer('mRNA_sPAcP',['eif4e'])
Monomer('S40', ['mrna'])
Monomer('S60', ['mrna'])

# 11) Protein products (sinks in your list)
Monomer('CycD')
Monomer('PSA')



# initial concentrations
Parameter('Lig_0', 100)
Parameter('Rec_0', 100)

Initial(EGF(r=None), Lig_0)
Initial(EGFR(l=None, d=None,grb2_shc=None, state='u'), Rec_0)

# rules
Parameter('kf_LR', 1)
Parameter('kr_LR', 1)
Rule('lig_binds_rec', EGF(r=None) + EGFR(l=None, d=None) | EGF(r=1) % EGFR(l=1, d=None), kf_LR, kr_LR)

Parameter('kf_Her2_cPAcP', 1)
Parameter('kr_Her2_cPAcP', 1)
Rule('Her2p_binds_cPAcP',Her2(d=ANY, cpacp=None,state='p') + cPAcP(her2=None) |
Her2(d=ANY, cpacp=1, state='p') % cPAcP(her2=1),
kf_Her2_cPAcP, kr_Her2_cPAcP)


# observable
Observable('Lig_free', EGF(r=None))
Observable('Lig_bound', EGF(r=ANY))
Observable('Rec_unphos', EGFR(state='u'))
Observable('Rec_phos', EGFR(state='p'))

if __name__ == '__main__':
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
