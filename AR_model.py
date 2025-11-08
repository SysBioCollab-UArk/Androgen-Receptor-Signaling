from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

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
Monomer('g_CycD', ['rnap','tf','ets_p','ap1_p'])
Monomer('g_PSA',  ['rnap','tf','ar_dimer'])

# 10) mRNAs and translation machinery
Monomer('eIF4E',     ['_4ebp1','mrna'])
Monomer('_40S',      ['mrna'])
Monomer('_60S',      ['_40s'])
Monomer('mRNA_cPAcP',['eif4e','_40s'])
Monomer('mRNA_sPAcP',['eif4e','_40s'])
Monomer('mRNA_CycD', ['eif4e','_40s'])
Monomer('mRNA_PSA',  ['eif4e','_40s'])


# 11) Protein products (sinks in your list)
Monomer('CycD')
Monomer('PSA')


# initial concentrations
Parameter('Lig_0', 100)
Parameter('Rec_0', 100)

Initial(EGF(r=None, loc='mem'), Lig_0)
Initial(EGFR(l=None, d=None, grb2_shc=None, state='u', loc='mem'), Rec_0)


# rules
def transcription_rule(protein, kf_kr_kcat, tfs=None, k_on_off=None):

    if len(np.array(kf_kr_kcat).shape) == 1:
        kf_kr_kcat = [kf_kr_kcat]

    if len(np.array(tfs).shape) == 0:
        tfs = [tfs]
        k_on_off = [k_on_off]

    gene = model.monomers['g_%s' % protein.name]
    mrna = model.monomers['mRNA_%s' % protein.name]

    kf, kr, kcat = \
        [Parameter('%s_%s_RNAp' % (k, gene.name), kf_kr_kcat[0][i]) for i, k in enumerate(['kf', 'kr', 'kcat'])]
    Rule('%s_binds_RNAp' % gene.name,
         gene(rnap=None, tf=None) + RNAp(gene=None) | gene(rnap=1, tf=None) % RNAp(gene=1), kf, kr)
    Rule('%s_RNAp_transcribe' % gene.name,
         gene(rnap=1, tf=None) % RNAp(gene=1) >> gene(rnap=1, tf=None) % RNAp(gene=1) + mrna(eif4e=None,_40s=None),
         kcat)

    if tfs is not None:
        for i, tf in enumerate(tfs):
            k_on, k_off = [Parameter('%s_%s_%s' % (k, gene.name, tf.monomer.name), k_on_off[i][j])
                           for j, k in enumerate(['kon', 'koff'])]
            Rule('%s_binds_%s' % (gene.name, tf.monomer.name),
                 gene(rnap=None, tf=None) + tf(gene=None) | gene(rnap=None, tf=1) % tf(gene=1), k_on, k_off)

            kf, kr, kcat = [Parameter('%s_%s_%s_RNAp' % (k, gene.name, tf.monomer.name), kf_kr_kcat[i+1][j])
                                      for j, k in enumerate(['kf', 'kr', 'kcat'])]
            Rule('%s_%s_binds_RNAp' % (gene.name, tf.monomer.name),
                 gene(rnap=None, tf=1) % tf(gene=1) + RNAp(gene=None) | gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2),
                 kf, kr)
            Rule('%s_%s_RNAp_transcribe' % (gene.name, tf.monomer.name),
                 gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2) >> gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2)
                 + mrna(eif4e=None, _40s=None), kcat)


transcription_rule(CycD, [[1, 1, 1], [1, 1, 1]], tfs=ETS(erk=None,gene=None,state='p'), k_on_off=[1, 1])

print(model.rules)

quit()

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
