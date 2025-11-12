from pysb import *
from pysb import MonomerPattern, as_complex_pattern
from util import set_model, create_transcription_rules, create_translation_rules
from itertools import product as cartesian_product

# Reimplementation of the androgen receptor signaling model from:
# Tasseff R, Nayak S, Salim S, Kaushik P, Rizvi N, Varner JD (2010)
# "Analysis of the Molecular Networks in Androgen Dependent and Independent Prostate Cancer
# Revealed Fragile and Robust Subsystems". PLoS ONE 5(1): e8864.
# https://doi.org/10.1371/journal.pone.0008864

Model()
set_model(model)  # put the model into the global namespace


# === MONOMERS ===

# 1) Ligands & Receptors
Monomer('EGF',  ['r','loc'], {'loc': ['extra','intra']})
Monomer('EGFR', ['l','d','grb2_shc','state','loc'],{'state':['u','p'], 'loc': ['extra','intra']})
Monomer('Her2', ['d','grb2_shc','cpacp','state'],{'state':['u','p']})
Monomer('Grb2', ['sos','shc','egfr_her2'])
Monomer('Sos',  ['grb2','ras_erk','pi3k'])
Monomer('Ras',  ['sos','gap','raf','state'],{'state':['GDP','GTP']})
Monomer('Shc',  ['r1', 'r2','grb2','state'],{'state':['u','p']})

# 2) Phosphatases acting on HER2/AR modules
Monomer('cPAcP',  ['r1', 'r2'])
Monomer('sPAcP',  ['r1', 'r2', 'loc'], {'loc': ['intra','extra']})

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
Monomer('g_sPAcP',['rnap','tf1', 'tf2'])
Monomer('g_CycD', ['rnap','tf'])
Monomer('g_PSA',  ['rnap','tf1', 'tf2'])

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

# Initial(EGF(r=None, loc='extra'), EGF_0)
# Initial(EGFR(l=None, d=None, grb2_shc=None, state='u', loc='extra'), EGFR_0)

# Monomer('sPAcP',  ['r1', 'r2', 'loc'], {'loc': ['intra','extra']})
# Monomer('T',      ['b','loc'], {'loc': ['intra', 'extra']})  # testosterone
# Monomer('mRNA_cPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
# Monomer('mRNA_sPAcP',['eif4e','_40s','elong'], {'elong': ['i','a']})
# Monomer('mRNA_CycD', ['eif4e','_40s','elong'], {'elong': ['i','a']})
# Monomer('mRNA_PSA',  ['eif4e','_40s','elong'], {'elong': ['i','a']})

init_params = {
    # Monomer('EGF',  ['r','loc'], {'loc': ['extra','intra']})
    'EGF_loc_extra_0': 100,  # 8 nM
    # Pase7 24.71 ± 23.66
    'Pase7_0': 24.71,
    # AR 192.40 ± 260.66
    # Monomer('AR',     ['lig','ar','erk','pase5','gene','state'], {'state':['u','p']})
    'AR_state_u_0': 192.40,
    # HSP 486.15 ± 659.55
    'HSP_0': 486.15,
    # Rase5a 81.11 ± 66.98
    'Rase5a_0': 81.11,
    # Her2 131.57 ± 111.41
    # Monomer('Her2', ['d','grb2_shc','cpacp','state'],{'state':['u','p']})
    'Her2_state_u_0': 131.57,
    # EGFR 115.41 ± 82.19
    # Monomer('EGFR', ['l','d','grb2_shc','state','loc'],{'state':['u','p'], 'loc': ['extra','intra']})
    'EGFR_state_u_loc_extra_0': 115.41,
    # Shc 91.19 ± 71.89
    # Monomer('Shc',  ['r1', 'r2','grb2','state'],{'state':['u','p']})
    'Shc_state_u_0': 91.19,
    # Grb2 80.32 ± 67.10
    'Grb2_0': 80.32,
    # Sos 35.79 ± 29.37
    'Sos_0': 35.79,
    # Ras-GDP 233.39 ± 456.75
    # Monomer('Ras',  ['sos','gap','raf','state'],{'state':['GDP','GTP']})
    'Ras_state_GDP_0': 233.39,
    # Raf 76.83 ± 54.31
    # Monomer('Raf',    ['ras','mek','pase1','state'],{'state':['u','p']})
    'Raf_state_u_0': 76.83,
    # MEK 1572.31 ± 2260.09
    # Monomer('MEK',    ['raf','erk','pase2','state'],{'state':['u','p','pp']})
    'MEK_state_u_0': 1572.31,
    # ERK 587.24 ± 401.24
    # Monomer('ERK',    ['ar', 'mek','pase3','sos','ets','ap1','state'],{'state':['u','p','pp']})
    'ERK_state_u_0': 587.24,
    # ETS 133.52 ± 150.03
    # Monomer('ETS',    ['erk','gene','state'],{'state':['u','p']})
    'ETS_state_u_0': 133.52,
    # AP1 107.34 ± 172.51
    # Monomer('AP1',    ['erk','gene','state'],{'state':['u','p']})
    'AP1_state_u_0': 107.34,
    # Pase1 181.67 ± 513.24
    'Pase1_0': 181.67,
    # Pase2 20.88 ± 11.78
    'Pase2_0': 20.88,
    # Pase3 22.76 ± 12.44
    'Pase3_0': 22.76,
    # Pase5 65.13 ± 82.13
    'Pase5_0': 65.13,
    # Pase6 168.28 ± 234.48
    'Pase6_0': 168.28,
    # GAP 60.71 ± 100.62
    'GAP_0': 60.71,
    # PI3K 174.54 ± 240.45
    # Monomer('PI3K',   ['egfr_her2','ptdins2','sos','state'],{'state':['i','act']})
    'PI3K_state_i_0': 174.54,
    # PtdIns2 131.27 ± 122.14
    'PtdIns2_0': 131.27,
    # PtdIns3 119.05 ± 94.11
    'PtdIns3_0': 119.05,
    # PTEN 123.84 ± 142.99
    'PTEN_0': 123.84,
    # Akt 332.36 ± 585.30
    # Monomer('Akt',    ['ptdins3','pdk1','tor','pase7','state'],{'state':['i','m','act']})
    'Akt_state_i_0': 332.36,
    # Pdk1 190.88 ± 237.37
    # Monomer('Pdk1',   ['ptdins3','akt','state'],{'state':['i','m']})
    'Pdk1_state_i_0': 190.88,
    # TOR 121.80 ± 120.53
    # Monomer('TOR',    ['akt','e4ebp1','state'],{'state':['i','act']})
    'TOR_state_i_0': 121.80,
    # 4E-BP1 136.67 ± 107.57
    # Monomer('_4EBP1', ['eif4e','tor','state'],{'state':['u','p']})
    '_4EBP1_state_u_0': 136.67,
    # eIF4E 3707.42 ± 3178.77
    'eIF4E_0': 3707.42,
    # g-PSA 3.29 ± 1.87
    'g_PSA_0': 3.29,
    # g-CycD 2.90 ± 4.98
    'g_CycD_0': 2.90,
    # g-cPAcP 0.09 ± 0.06
    'g_cPAcP_0': 0.09,
    # g-sPAcP 0.11 ± 0.09
    'g_sPAcP_0': 0.11,
    # RNAp 371.62 ± 314.59
    'RNAp_0': 371.62,
    # 40S 8203.19 ± 16956.30
    '_40S_0': 8203.19,
    # 60S 4732.73 ± 4700.66
    '_60S_0': 4732.73
}

# create initials for all monomers and all states
for mon in model.monomers:
    sites_NONE = [(site, None) for site in mon.sites if site not in mon.site_states.keys()]
    for states in list(cartesian_product(*[list(mon.site_states[site]) for site in mon.site_states.keys()])):
        sites_STATES = []
        suffix = ''
        for i, site in enumerate(mon.site_states.keys()):
            sites_STATES.append((site, states[i]))
            suffix += '_%s_%s' % (site, states[i])
        mp = MonomerPattern(mon, dict(sites_NONE + sites_STATES), None)
        pname = '%s%s_0' % (mon.name, suffix)
        model.add_initial(Initial(as_complex_pattern(mp), Parameter(pname, init_params.get(pname, 0))))


# === RULES ===

# EGF+EGFR	↔	EGFR-EGF
Parameter('kf_EGF_binds_EGFR', 2.215)
Parameter('kr_EGF_binds_EGFR', 1.343e-3)
Rule('EGF_binds_EGFR',
     EGF(r=None, loc='extra') + EGFR(l=None, d=None, grb2_shc=None, state='u', loc='extra') |
     EGF(r=1, loc='extra') % EGFR(l=1, d=None, grb2_shc=None, state='u', loc='extra'),
     kf_EGF_binds_EGFR, kr_EGF_binds_EGFR)

# 2*EGFR-EGF	↔	EGFR-EGF-2
Parameter('kf_EGFR_EGF_dimerize', 0.3701)
Parameter('kr_EGFR_EGF_dimerize', 0.1708)
Rule('EGFR_EGF_dimerization',
     EGF(r=1, loc='extra') % EGFR(l=1, d=None, grb2_shc=None, state='u', loc='extra') +
     EGF(r=2, loc='extra') % EGFR(l=2, d=None, grb2_shc=None, state='u', loc='extra') |
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=None, state='u', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=None, state='u', loc='extra'),
     kf_EGFR_EGF_dimerize, kr_EGFR_EGF_dimerize)

# EGFR-EGF-2	↔	EGFR-EGF-2-p
Parameter('kf_EGFR_EGF_phos', 1.864)
Parameter('kr_EGFR_EGF_phos', 2.707e-2)
Rule('EGFR_autophosphorylation',
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=None, state='u', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=None, state='u', loc='extra') |
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=None, state='p', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=None, state='p', loc='extra'),
     kf_EGFR_EGF_phos, kr_EGFR_EGF_phos)

# 2*Her2	↔	Her2-2
Parameter('kf_Her2_dimer', 4.98e-2)
Parameter('kr_Her2_dimer', 0.1756)
Rule('Her2_dimerization',
     Her2(d=None, grb2_shc=None, cpacp=None, state='u') + Her2(d=None, grb2_shc=None, cpacp=None, state='u') |
     Her2(d=1, grb2_shc=None, cpacp=None, state='u') % Her2(d=1, grb2_shc=None, cpacp=None, state='u'),
     kf_Her2_dimer, kr_Her2_dimer)

# Her2-2	↔	Her2-2-p
Parameter('kf_Her2_2_phos', 2.032e-2)
Parameter('kr_Her2_2_phos', 1.472e-5)
Rule('Her2_2_phosphorylation',
     Her2(d=1, grb2_shc=None, cpacp=None, state='u') % Her2(d=1, grb2_shc=None, cpacp=None, state='u') |
     Her2(d=1, grb2_shc=None, cpacp=None, state='p') % Her2(d=1, grb2_shc=None, cpacp=None, state='p'),
     kf_Her2_2_phos, kr_Her2_2_phos)

# EGFR-EGF-2-p+Grb2	↔	EGFR-EGF-2-p-Grb2	1.068E0±3.282E0	7.018E-1±4.517E-1	-
# TODO ...

# EGFR-EGF-2-p-Grb2+Sos	↔	EGFR-EGF-2-p-Grb2-Sos	4.244E-1±6.357E-1	3.506E0±5.793E0	-
# TODO ...

# EGFR-EGF-2-p-Grb2-Sos	↔	EGFR-EGF-2-p+Grb2-Sos	1.159E0±2.35E0	5.034E-2±6.075E-2	-
# TODO ...

# Grb2+Sos	↔	Grb2-Sos	4.104E-3, 2.548E-4
Parameter('kf_Grb2_binds_Sos', 4.104E-3)
Parameter('kr_Grb2_binds_Sos', 2.548E-4)
Rule('Grb2_binds_Sos',
     Grb2(sos=None, shc=None, egfr_her2=None) + Sos(grb2=None, ras_erk=None, pi3k=None) |
     Grb2(sos=1, shc=None, egfr_her2=None) % Sos(grb2=1, ras_erk=None, pi3k=None),
     kf_Grb2_binds_Sos, kr_Grb2_binds_Sos)

# EGFR-EGF-2-p-Grb2-Sos+Ras-GDP	↔	EGFR-EGF-2-p-Grb2-Sos-Ras-GDP	2.183E-2±1.842E-2	8.377E-1±9.551E-1	-
# EGFR-EGF-2-p-Grb2-Sos-Ras-GDP	→	EGFR-EGF-2-p-Grb2-Sos+Ras-GTP	-	-	4.28E0±3.525E0
# TODO ...

# Her2-2-p+Grb2	↔	Her2-2-p-Grb2	1.976E-2±2.734E-2	9.558E-1±8.677E-1	-
# TODO ...

# Her2-2-p-Grb2+Sos	↔	Her2-2-p-Grb2-Sos	2.395E-1±2.953E-1	1.49E0±2.296E0	-
# TODO ...

# Her2-2-p-Grb2-Sos	↔	Her2-2-p+Grb2-Sos	3.623E-2±2.45E-2	5.343E-2±5.586E-2	-
# TODO ...

# Her2-2-p-Grb2-Sos+Ras-GDP	↔	Her2-2-p-Grb2-Sos-Ras-GDP	3.053E-2±1.959E-2	4.568E-1±3.063E-1	-
# Her2-2-p-Grb2-Sos-Ras-GDP	→	Her2-2-p-Grb2-Sos+Ras-GTP	-	-	5.685E0±5.969E0
# TODO ...

# EGFR-EGF-2-p+Shc	↔	EGFR-EGF-2-p-Shc
# EGFR-EGF-2-p-Shc	→	EGFR-EGF-2-p-Shc-p
Parameter('kf_EGFR_EFG_2_p_binds_Shc', 0.6682)
Parameter('kr_EGFR_EFG_2_p_binds_Shc', 2.79)
Parameter('kcat_EGFR_EGF_2_p_binds_Shc', 13.27)
Rule('EGFR_EGF_2_p_binds_Shc',
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=None, state='p', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=None, state='p', loc='extra') +
     Shc(r1=None, r2=None, grb2=None, state='u') |
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=4, state='p', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=5, state='p', loc='extra') %
     Shc(r1=4, r2=5, grb2=None, state='u'),
     kf_EGFR_EFG_2_p_binds_Shc, kr_EGFR_EFG_2_p_binds_Shc)
Rule('EGFR_EGF_2_p_Shc_phos_Shc_p',
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=4, state='p', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=5, state='p', loc='extra') %
     Shc(r1=4, r2=5, grb2=None, state='u') >>
     EGF(r=1, loc='extra') % EGFR(l=1, d=3, grb2_shc=4, state='p', loc='extra') %
     EGF(r=2, loc='extra') % EGFR(l=2, d=3, grb2_shc=5, state='p', loc='extra') %
     Shc(r1=4, r2=5, grb2=None, state='p'),
     kcat_EGFR_EGF_2_p_binds_Shc)

# EGFR-EGF-2-p-Shc-p	↔	EGFR-EGF-2-p+Shc-p	1.464E0±1.058E0	6.464E-3±8.466E-3	-
# TODO ...

# Shc-p	→	Shc	-	-	5.153E0±5.371E0
# TODO ...

# EGFR-EGF-2-p-Shc-p+Grb2	↔	EGFR-EGF-2-p-Shc-p-Grb2	6.308E-2±4.851E-2	8.773E-1±6.182E-1	-
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2+Sos	↔	EGFR-EGF-2-p-Shc-p-Grb2-Sos	2.065E-1±3.931E-1	2.176E-1±2.282E-1	-
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2-Sos	↔	EGFR-EGF-2-p+Shc-p-Grb2-Sos	9.887E-1±1.048E0	1.562E-3±1.142E-3	-
# TODO ...

# Shc-p+Grb2-Sos	↔	Shc-p-Grb2-Sos
Parameter('kf_Shc_p_binds_Grb2_Sos', 3.755e-2)
Parameter('kr_Shc_p_binds_Grb2_Sos', 0.1613)
Rule('Shcp_binds_Grb2',
     Shc(r1=None, r2=None, grb2=None, state='p') +
     Grb2(sos=1, shc=None, egfr_her2=None) % Sos(grb2=1, ras_erk=None, pi3k=None) |
     Shc(r1=None, r2=None, grb2=2, state='p') %
     Grb2(sos=1, shc=2, egfr_her2=None) % Sos(grb2=1, ras_erk=None, pi3k=None),
     kf_Shc_p_binds_Grb2_Sos, kr_Shc_p_binds_Grb2_Sos)

# EGFR-EGF-2-p-Shc-p-Grb2-Sos+Ras-GDP	↔	EGFR-EGF-2-p-Shc-p-Grb2-Sos-Ras-GDP	3.412E-2±3.79E-2	2.066E-1±2.164E-1	-
# EGFR-EGF-2-p-Shc-p-Grb2-Sos-Ras-GDP	→	EGFR-EGF-2-p-Shc-p-Grb2-Sos+Ras-GTP	-	-	4.399E0±7.045E0
# TODO ...

# Her2-2-p+Shc	↔	Her2-2-p-Shc	1.208E-1±1.537E-1	9.921E0±2.087E1	-
# Her2-2-p-Shc	→	Her2-2-p-Shc-p	-	-	2.995E1±6.778E1
# TODO ...

# Her2-2-p-Shc-p	↔	Her2-2-p+Shc-p	4.703E0±1.342E1	3.983E-3±4.616E-3	-
# TODO ...

# Her2-2-p-Shc-p+Grb2	↔	Her2-2-p-Shc-p-Grb2	5.805E-2±7.916E-2	1.127E0±1.51E0	-
# TODO ...

# Her2-2-p-Shc-p-Grb2+Sos	↔	Her2-2-p-Shc-p-Grb2-Sos	1.359E-1±1.23E-1	2.241E0±7.899E0	-
# TODO ...

# Her2-2-p-Shc-p-Grb2-Sos	↔	Her2-2-p+Shc-p-Grb2-Sos	9.861E-1±1.159E0	1.294E-2±3.787E-2	-
# TODO ...

# Her2-2-p-Shc-p-Grb2-Sos+Ras-GDP	↔	Her2-2-p-Shc-p-Grb2-Sos-Ras-GDP	1.749E-2±1.219E-2	6.129E-1±1.44E0	-
# Her2-2-p-Shc-p-Grb2-Sos-Ras-GDP	→	Her2-2-p-Shc-p-Grb2-Sos+Ras-GTP	-	-	1.882E1±4.063E1
# TODO ...

# Her2-2-p+cPAcP	↔	Her2-2-p-cPAcP
# Her2-2-p-cPAcP	→	Her2-2+cPAcP
Parameter('kf_Her2_2_p_cPAcP_dephos', 17.07)
Parameter('kr_Her2_2_p_cPAcP_dephos', 5.325e-2)
Parameter('kcat_Her2_2_p_cPAcP_dephos', 20.12)
Rule('Her2_2_p_binds_cPAcP',
     Her2(d=1, grb2_shc=None, cpacp=None, state='p') % Her2(d=1, grb2_shc=None, cpacp=None, state='p') +
     cPAcP(r1=None, r2=None) |
     Her2(d=1, grb2_shc=None, cpacp=2, state='p') % Her2(d=1, grb2_shc=None, cpacp=3, state='p') %
     cPAcP(r1=2, r2=3),
     kf_Her2_2_p_cPAcP_dephos, kr_Her2_2_p_cPAcP_dephos)
Rule('Her2_2_p_cPAcP_dephos',
     Her2(d=1, grb2_shc=None, cpacp=2, state='p') % Her2(d=1, grb2_shc=None, cpacp=3, state='p') %
     cPAcP(r1=2, r2=3) >>
     Her2(d=1, grb2_shc=None, cpacp=None, state='u') % Her2(d=1, grb2_shc=None, cpacp=None, state='u') +
     cPAcP(r1=None, r2=None),
     kcat_Her2_2_p_cPAcP_dephos)

# Her2-2+sPAcP	↔	Her2-2-sPAcP	8.951E0±9.414E0	1.248E-3±7.101E-4	-
# Her2-2-sPAcP	→	Her2-2-p+sPAcP	-	-	2.288E1±2.252E1
# TODO ...

# Ras-GTP+GAP	↔	Ras-GTP-GAP	1.032E-1±1.526E-1	1.149E0±1.23E0	-
# Ras-GTP-GAP	→	Ras-GDP+GAP	-	-	4.785E-1±3.57E-1
# TODO ...

# Ras-GTP+Raf	↔	Ras-GTP-Raf	5.455E-3±4.49E-3	2.097E-2±1.256E-2	-
# Ras-GTP-Raf	→	Ras-GTP+Raf-p	-	-	5.858E0±6.133E0
# TODO ...

# Raf-p+Pase1	↔	Raf-p-Pase1	5.166E-1±5.821E-1	1.77E0±1.446E0	-
# Raf-p-Pase1	→	Raf+Pase1	-	-	3.978E0±2.632E0
# TODO ...

# MEK+Raf-p	↔	MEK-Raf-p	6.055E-2±5.209E-2	9.006E-2±1.189E-1	-
# MEK-Raf-p	→	MEK-p+Raf-p	-	-	1.409E1±2.857E1
# TODO ...

# MEK-p+Raf-p	↔	MEK-p-Raf-p	2.145E-1±6.272E-1	1.056E-1±1.282E-1	-
# MEK-p-Raf-p	→	MEK-pp+Raf-p	-	-	2.949E0±2.16E0
# TODO ...

# ERK+MEK-pp	↔	ERK-MEK-pp	1.677E-3±1.72E-3	4.956E-1±3.873E-1	-
# ERK-MEK-pp	→	ERK-p+MEK-pp	-	-	1.095E1±1.112E1
# TODO ...

# ERK-p+MEK-pp	↔	ERK-p-MEK-pp	2.09E-3±1.47E-3	3.124E0±8.012E0	-
# ERK-p-MEK-pp	→	ERK-pp+MEK-pp	-	-	8.435E0±6.538E0
# TODO ...

# MEK-p+Pase2	↔	MEK-p-Pase2	1.04E-3±1.072E-3	7.569E0±1.607E1	-
# MEK-p-Pase2	→	MEK+Pase2	-	-	7.221E-1±8.541E-1
# TODO ...

# MEK-pp+Pase2	↔	MEK-pp-Pase2	6.319E-2±4.223E-2	4.293E0±2.501E0	-
# MEK-pp-Pase2	→	MEK-p+Pase2	-	-	2.617E-2±2.29E-2
# TODO ...

# ERK-p+Pase3	↔	ERK-p-Pase3	2.408E0±2.499E0	3.366E-1±2.34E-1	-
# ERK-p-Pase3	→	ERK+Pase3	-	-	2.7E0±4.156E0
# TODO ...

# ERK-pp+Pase3	↔	ERK-pp-Pase3	5.008E-2±3.168E-2	5.499E0±8.337E0	-
# ERK-pp-Pase3	→	ERK-p+Pase3	-	-	2.789E0±5.206E0
# TODO ...

# EGFR-EGF-2-p-Grb2-Sos+ERK-pp	↔	EGFR-EGF-2-p-Grb2-Sos-ERK-pp	2.239E0±2.85E0	6.502E-4±5.11E-4	-
# EGFR-EGF-2-p-Grb2-Sos-ERK-pp	→	EGFR-EGF-2-p-Grb2+Sos+ERK-pp	-	-	1.351E0±1.222E0
# TODO ...

# Her2-2-p-Grb2-Sos+ERK-pp	↔	Her2-2-p-Grb2-Sos-ERK-pp	1.856E0±1.211E0	1.559E-1±2.128E-1	-
# Her2-2-p-Grb2-Sos-ERK-pp	→	Her2-2-p-Grb2+Sos+ERK-pp	-	-	3.075E0±3.738E0
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2-Sos+ERK-pp	↔	EGFR-EGF-2-p-Shc-p-Grb2-Sos-ERK-pp	1.185E0±1.099E0	8.315E-4±6.405E-4	-
# EGFR-EGF-2-p-Shc-p-Grb2-Sos-ERK-pp	→	EGFR-EGF-2-p-Shc-p-Grb2+Sos+ERK-pp	-	-	1.231E1±4.931E1
# TODO ...

# Her2-2-p-Shc-p-Grb2-Sos+ERK-pp	↔	Her2-2-p-Shc-p-Grb2-Sos-ERK-pp	3.378E0±5.412E0	3.747E-1±8.356E-1	-
# Her2-2-p-Shc-p-Grb2-Sos-ERK-pp	→	Her2-2-p-Shc-p-Grb2+Sos+ERK-pp	-	-	5.262E0±9.916E0
# TODO ...

# EGF	→	EGFi	-	-	0E0±0E0
# TODO ...

# EGFR	↔	EGFRi	1.179E-2±1.056E-2	1.599E-1±2.308E-1	-
# TODO ...

# EGFR-EGF	→	EGFR-EGFi	-	-	1.579E-1±3.334E-1
# TODO ...

# EGFRi+EGFi	↔	EGFR-EGFi	1.041E0±1.185E0	1.251E-1±4.785E-1	-
# TODO ...

# EGFR-EGF-2	→	EGFR-EGF-2i	-	-	6.648E-3±4.899E-3
# TODO ...

# 2*EGFR-EGFi	↔	EGFR-EGF-2i	1.432E-2±1.301E-2	4.332E0±7.81E0	-
# TODO ...

# EGFR-EGF-2-p	→	EGFR-EGF-2-pi	-	-	7.456E-2±8.563E-2
# TODO ...

# EGFR-EGF-2i	↔	EGFR-EGF-2-pi	2.347E0±3.074E0	1.54E-1±1.493E-1	-
# TODO ...

# EGFR-EGF-2-p-Grb2	→	EGFR-EGF-2-pi-Grb2	-	-	7.874E-2±9.134E-2
# TODO ...

# EGFR-EGF-2-pi+Grb2	↔	EGFR-EGF-2-pi-Grb2	3.449E-3±2.48E-3	3.572E-1±3.045E-1	-
# TODO ...

# EGFR-EGF-2-p-Grb2-Sos	→	EGFR-EGF-2-pi-Grb2-Sos	-	-	1.541E-1±2.247E-1
# TODO ...

# EGFR-EGF-2-pi+Grb2-Sos	↔	EGFR-EGF-2-pi-Grb2-Sos	3.954E-3±3.183E-3	3.359E-1±2.71E-1	-
# TODO ...

# EGFR-EGF-2-pi-Grb2+Sos	↔	EGFR-EGF-2-pi-Grb2-Sos	1.13E-2±1.034E-2	6.348E-1±6.501E-1	-
# TODO ...

# EGFR-EGF-2-pi-Grb2-Sos+Ras-GDP	↔	EGFR-EGF-2-pi-Grb2-Sos-Ras-GDP	3.459E-3±3.064E-3	7.339E-2±8.78E-2	-
# EGFR-EGF-2-pi-Grb2-Sos-Ras-GDP	→	EGFR-EGF-2-pi-Grb2-Sos+Ras-GTP	-	-	1.199E0±9.562E-1
# TODO ...

# EGFR-EGF-2-pi-Grb2-Sos+ERK-pp	↔	EGFR-EGF-2-pi-Grb2-Sos-ERK-pp	1.427E0±1.503E0	3.947E-4±2.363E-4	-
# EGFR-EGF-2-pi-Grb2-Sos-ERK-pp	→	EGFR-EGF-2-pi-Grb2+Sos+ERK-pp	-	-	1.617E0±1.465E0
# TODO ...

# EGFR-EGF-2-pi+Shc	↔	EGFR-EGF-2-pi-Shc	1 .367E-1±1.947E-1	6.82E0±6.97E0	-
# TODO ...

# EGFR-EGF-2-p-Shc	→	EGFR-EGF-2-pi-Shc	-	-	1.652E-1±3.045E-1
# TODO ...

# EGFR-EGF-2-pi-Shc	→	EGFR-EGF-2-pi-Shc-p	-	-	9.176E0±1.724E1
# TODO ...

# EGFR-EGF-2-p-Shc-p	→	EGFR-EGF-2-pi-Shc-p	-	-	9.441E-2±6.697E-2
# TODO ...

# EGFR-EGF-2-pi-Shc-p	↔	EGFR-EGF-2-pi+Shc-p	6.97E0±1.16E1	1.036E-3±7.374E-4	-
# TODO ...

# EGFR-EGF-2-pi-Shc-p+Grb2	↔	EGFR-EGF-2-pi-Shc-p-Grb2	4.201E-3±7.391E-3	2.866E0±6.449E0	-
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2	→	EGFR-EGF-2-pi-Shc-p-Grb2	-	-	1.38E-1±1.346E-1
# TODO ...

# EGFR-EGF-2-pi-Shc-p-Grb2+Sos	↔	EGFR-EGF-2-pi-Shc-p-Grb2-Sos	9.048E-3±4.866E-3	5.588E-1±7.609E-1	-
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2-Sos	→	EGFR-EGF-2-pi-Shc-p-Grb2-Sos	-	-	4.926E-1±1.139E0
# TODO ...

# EGFR-EGF-2-pi-Shc-p-Grb2-Sos	↔	EGFR-EGF-2-pi+Shc-p-Grb2-Sos	3.537E0±3.537E0	3.606E-4±3.165E-4	-
# TODO ...

# EGFR-EGF-2-pi-Shc-p-Grb2-Sos+Ras-GDP	↔	EGFR-EGF-2-pi-Shc-p-Grb2-Sos-Ras-GDP	5.62E-3±4.616E-3	7.134E-1±1.236E0	-
# EGFR-EGF-2-pi-Shc-p-Grb2-Sos-Ras-GDP	→	EGFR-EGF-2-pi-Shc-p-Grb2-Sos+Ras-GTP	-	-	1.142E0±8.919E-1
# TODO ...

# EGFR-EGF-2-pi-Shc-p-Grb2-Sos+ERK-pp	↔	EGFR-EGF-2-pi-Shc-p-Grb2-Sos-ERK-pp	1.924E0±2.013E0	5.668E-4±4.798E-4	-
# EGFR-EGF-2-pi-Shc-p-Grb2-Sos-ERK-pp	→	EGFR-EGF-2-pi-Shc-p-Grb2+Sos+ERK-pp	-	-	4.336E0±4.556E0
# TODO ...

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
Parameter('kcat_AR_binds_T', 5.458e-2)
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
     AR(lig=None, ar=2, erk=None, pase5=None, gene=None, state='p') %
     AR(lig=1, ar=2, erk=None, pase5=None, gene=None, state='p') % T(b=1),
     kon_AR_p_AR_p_T, koff_AR_p_AR_p_T)

# 2*AR-p-T	↔	AR-p-T-2
Parameter('kon_AR_p_T_dimerize', 0.5835)
Parameter('koff_AR_p_T_dimerize', 5.149e-4)
Rule('AR_p_T_dimerizes',
     AR(lig=1, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=1) +
     AR(lig=2, ar=None, erk=None, pase5=None, gene=None, state='p') % T(b=2) |
     AR(lig=1, ar=3, erk=None, pase5=None, gene=None, state='p') % T(b=1) %
     AR(lig=2, ar=3, erk=None, pase5=None, gene=None, state='p') % T(b=2),
     kon_AR_p_T_dimerize, koff_AR_p_T_dimerize)

# 2*AR-p	↔	AR-p-2
Parameter('kon_AR_p_dimerize', 0.2848)
Parameter('koff_AR_p_dimerize', 0.1281)
Rule('AR_p_dimerizes',
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='p') +
     AR(lig=None, ar=None, erk=None, pase5=None, gene=None, state='p') |
     AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p') %
     AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'),
     kon_AR_p_dimerize, koff_AR_p_dimerize)

# AR+DHT	↔	AR-DHT	2.486E0±2.641E0	5.42E-5±6.636E-5	-
# AR-DHT	→	AR-p-DHT	5.436E-1±4.693E-1	-	-
# TODO ...

# AR-p-DHT+AR-p-T	↔	AR-p-DHT-AR-p-T	6.895E-1±1.445E0	7.915E-4±5.504E-4	-
# TODO ...

# AR-p-DHT+AR-p	↔	AR-p-DHT-AR-p	6.064E-1±1.283E0	3.362E-3±6.635E-3	-
# TODO ...

# 2*AR-p-DHT	↔	AR-p-DHT-2	1.026E0±1.066E0	1.013E-3±1.461E-3	-
# TODO ...

# AR-p+Pase5	↔	AR-p-Pase5	5.409E-4±5.152E-4	4.6E-3±3.726E-3	-
# AR-p-Pase5	→	AR+Pase5	-	-	3.637E-3±2.758E-3
# TODO ...

# AR-p-T+Pase5	↔	AR-p-T-Pase5	1.671E-3±2.186E-3	7.194E-3±4.179E-3	-
# AR-p-T-Pase5	→	AR-T+Pase5	-	-	5.853E-3±5.596E-3
# TODO ...

# AR-p-DHT+Pase5	↔	AR-p-DHT-Pase5	7.907E-4±1.029E-3	7.667E-3±7.267E-3	-
# AR-p-DHT-Pase5	→	AR-DHT+Pase5	-	-	5.328E-2±1.285E-1
# TODO ...

# ERK-pp+ETS	↔	ERK-pp-ETS	2.109E-3±3.444E-3	4.624E-1±3.654E-1	-
# ERK-pp-ETS	→	ERK-pp+ETS-p	-	-	2.534E-2±1.687E-2
# TODO ...

# ETS-p+Pase5	↔	ETS-p-Pase5	3.753E0±3.797E0	2.548E-3±6.034E-3	-
# ETS-p-Pase5	→	ETS+Pase5	-	-	8.124E0±9.856E0
# TODO ...

# ERK-pp+AP1	↔	ERK-pp-AP1	1.403E-3±1.096E-3	5.971E-1±4.652E-1	-
# ERK-pp-AP1	→	ERK-pp+AP1-p	-	-	2.556E-2±2.661E-2
# TODO ...

# AP1-p+Pase6	↔	AP1-p-Pase6	8.022E0±1.209E1	8.007E-4±7.727E-4	-
# AP1-p-Pase6	→	AP1+Pase6	-	-	1.54E1±2.376E1
# TODO ...

# Her2-2-p-Grb2-Sos+PI3K	↔	Her2-2-p-Grb2-Sos-PI3K	2.125E-1±2.91E-1	1.412E-2±3.191E-2	-
# Her2-2-p-Grb2-Sos-PI3K	→	Her2-2-p-Grb2-Sos+Act-PI3K	-	-	1.941E-1±3.754E-1
# TODO ...

# Her2-2-p-Shc-p-Grb2-Sos+PI3K	↔	Her2-2-p-Shc-p-Grb2-Sos-PI3K	8.255E-2±8.716E-2	5.049E-3±7.031E-3	-
# Her2-2-p-Shc-p-Grb2-Sos-PI3K	→	Her2-2-p-Shc-p-Grb2-Sos+Act-PI3K	-	-	1.93E-1±1.462E-1
# TODO ...

# EGFR-EGF-2-p-Grb2-Sos+PI3K	↔	EGFR-EGF-2-p-Grb2-Sos-PI3K	4.34E-1±1.47E0	1.849E-2±4.73E-2	-
# EGFR-EGF-2-p-Grb2-Sos-PI3K	→	EGFR-EGF-2-p-Grb2-Sos+Act-PI3K	-	-	1.684E-1±2.04E-1
# TODO ...

# EGFR-EGF-2-p-Shc-p-Grb2-Sos+PI3K	↔	EGFR-EGF-2-p-Shc-p-Grb2-Sos-PI3K	1.722E-1±1.503E-1	8.976E-3±1.111E-2	-
# EGFR-EGF-2-p-Shc-p-Grb2-Sos-PI3K	→	EGFR-EGF-2-p-Shc-p-Grb2-Sos+Act-PI3K	-	-	1.01E-1±1.081E-1
# TODO ...

# PtdIns2+Act-PI3K	↔	PtdIns2-Act-PI3K	1.983E-1±1.959E-1	1.56E-2±9.585E-3	-
# PtdIns2-Act-PI3K	→	PtdIns3+Act-PI3K	-	-	9.81E-2±5.877E-2
# TODO ...

# PtdIns3+PTEN	↔	PtdIns3-PTEN	3.036E-1±3.942E-1	2.262E-2±3.503E-2	-
# PtdIns3-PTEN	→	PtdIns2+PTEN	-	-	3.081E-1±2.74E-1
# TODO ...

# PtdIns3+Akt	↔	PtdIns3-Akt	4.11E-1±7.091E-1	7.744E-3±1.008E-2	-
# PtdIns3-Akt	→	PtdIns3+Akt-m	-	-	2.839E-1±3.552E-1
# TODO ...

# PtdIns3+Pdk1	↔	PtdIns3-Pdk1	2.436E-1±4.38E-1	4.369E-3±7.769E-3	-
# PtdIns3-Pdk1	→	PtdIns3+Pdk1-m	-	-	7.83E0±2.249E1
# TODO ...

# Pdk1-m+Akt-m	↔	Pdk1-m-Akt-m	9.893E-2±6.529E-2	1.874E-2±4.907E-2	-
# Pdk1-m-Akt-m	→	Pdk1+Act-Akt	-	-	2.096E-1±3.052E-1
# TODO ...

# Act-Akt+TOR	↔	Act-Akt-TOR	1.389E-1±1.405E-1	1.102E-1±7.477E-2	-
# Act-Akt-TOR	→	Akt+Act-TOR	-	-	2.551E-1±1.76E-1
# TODO ...

# 4E-BP1+eIF4E	↔	4E-BP1-eIF4E
Parameter('kf_4EBP1_binds_eIF4E', 0.1779)
Parameter('kr_4EBP1_binds_eIF4E', 0.158)
Rule('_4EBP1_binds_eIF4E',
     _4EBP1(eif4e=None, tor=None, state='u') + eIF4E(mrna_4ebp1=None) |
     _4EBP1(eif4e=1, tor=None, state='u') % eIF4E(mrna_4ebp1=1),
     kf_4EBP1_binds_eIF4E, kr_4EBP1_binds_eIF4E)

# 4E-BP1+Act-TOR	↔	4E-BP1-Act-TOR	1.347E-1±1.643E-1	1.947E-1±2.04E-1	-
# 4E-BP1-Act-TOR	→	4E-BP1-P+Act-TOR	-	-	2.5E-1±2.358E-1
# TODO ...

# T-e	↔	T
Parameter('kf_T_extra_to_intra', 1.449)
Parameter('kr_T_extra_to_intra', 1.555)
Rule('T_extra_to_intra',
     T(b=None, loc='extra') | T(b=None, loc='intra'), kf_T_extra_to_intra, kr_T_extra_to_intra)

# sPAcP	→	sPAcP-e
Parameter('k_sPAcP_intra_to_extra', 1.86)
Rule('spAcP_intra_to_extra',
     sPAcP(r1=None, r2=None, loc='intra') >> sPAcP(r1=None, r2=None, loc='extra'), k_sPAcP_intra_to_extra)

# EGFi	→	[]	-	-	1.057E0±1.053E0
# TODO ...

# 2*cPAcP	↔	cPAcP-2	8.195E-2±1.868E-1	9.026E-2±6.808E-2	-
# TODO ...

# 2*cPAcP-2	↔	cPAcP-4	8.039E-2±1.344E-1	6.186E-2±3.229E-2	-
# TODO ...

# 2*Her2-2-p+cPAcP-2	↔	2Her2-2-p-cPAcP-2	1.292E2±3.518E2	2.069E-1±2.448E-1	-
# 2Her2-2-p-cPAcP-2	→	2*Her2-2+cPAcP-2	-	-	9.256E0±6.196E0
# TODO ...

# 4*Her2-2-p+cPAcP-4	↔	4Her2-2-p-cPAcP-4	1.306E1±1.071E1	1.127E-2±1.019E-2	-
# 4Her2-2-p-cPAcP-4	→	4*Her2-2+cPAcP-4	-	-	7.811E0±5.607E0
# TODO ...

# Act-Akt+Pase7	↔	Act-Akt-Pase7	1.765E-3±1.65E-3	1.226E-3±2.158E-3	-
# Act-Akt-Pase7	→	Akt+Pase7	-	-	1.861E-3±3.179E-3
# TODO ...

# 4E-BP1-eIF4E+Act-TOR	↔	4E-BP1-eIF4E-Act-TOR	1.603E-3±1.851E-3	7.501E-4±5.262E-4	-
# 4E-BP1-eIF4E-Act-TOR	→	4E-BP1-P+eIF4E+Act-TOR	-	-	1.933E-3±4.962E-3
# TODO ...


# === TRANSCRIPTION/TRANSLATION RULES ===

# *** cPAcP transcription ***
# g-cPAcP+RNAp	↔	g-cPAcP-RNAp
# g-cPAcP-RNAp	→	g-cPAcP+RNAp+mRNA-cPAcP
kf_kr_kcat = [[5.788e-2, 7.556e-2, 1.311e-2]]  # basal
# mRNA-cPAcP	→	[]
k_mRNA_deg = 1.3
# AR-p-2+g-cPAcP	↔	AR-p-2-g-cPAcP	1.618E-3±1.482E-3	9.872E-5±6.594E-5
# AR-p-DHT-2+g-cPAcP	↔	AR-p-DHT-2-g-cPAcP	8.423E-4±1.136E-3	1.366E-7±1.582E-7
# AR-p-DHT-AR-p+g-cPAcP	↔	AR-p-DHT-AR-p-g-cPAcP	4.544E-3±3.674E-3	1.314E-7±1.111E-7
# AR-p-T-2+g-cPAcP	↔	AR-p-T-2-g-cPAcP	6.374E-2±3.979E-2	2.954E-4±4.596E-4
# AR-p-AR-p-T+g-cPAcP	↔	AR-p-AR-p-T-g-cPAcP	6.191E-2±4.34E-2	1.099E-4±1.07E-4
# AR-p-DHT-AR-p-T+g-cPAcP	↔	AR-p-DHT-AR-p-T-g-cPAcP	1.543E-1±1.874E-1	2.343E-4±3.372E-4
tf_names_AR = ['AR_p_2', 'AR_p_DHT_2', 'AR_p_DHT_AR_p', 'AR_p_T_2', 'AR_p_T_AR_p', 'AR_p_DHT_AR_p_T']
tf_species_AR = [
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'),
    #
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    DHT(b=3) % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p'),
    #
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'),
    #
    T(b=2, loc='intra') % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    T(b=3, loc='intra') % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p'),
    #
    T(b=2, loc='intra') % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    AR(lig=None, ar=1, erk=None, pase5=None, gene=None, state='p'),
    #
    DHT(b=2) % AR(lig=2, ar=1, erk=None, pase5=None, gene=None, state='p') %
    T(b=3, loc='intra') % AR(lig=3, ar=1, erk=None, pase5=None, gene=None, state='p')
]
k_tf_on_off = [[1.618e-3, 9.872e-5], [8.423e-4, 1.366e-7], [4.544e-3, 1.314e-7], [6.374e-2, 2.954e-4],
               [6.191e-2, 1.099e-4], [0.1543,   2.343e-4]]
create_transcription_rules(cPAcP, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species_AR, tf_names_AR)),
                           k_tf_on_off=k_tf_on_off)

# *** cPAcP translation ***
# mRNA-cPAcP+eIF4E	↔	mRNA-cPAcP-eIF4E
# mRNA-cPAcP-eIF4E+40S	↔	mRNA-cPAcP-eIF4E-40S
# mRNA-cPAcP-eIF4E-40S+60S	↔	mRNA-cPAcP-eIF4E-40S-60S
kf_kr = [
    [1.351e-2, 9.892e-3],  # eIF4E binding
    [0.233, 9.908e-3],     # 40S binding
    [0.9053, 1.499e-3]     # 60S binding
]
# mRNA-cPAcP-eIF4E-40S-60S	→	Rm-cPAcP+eIF4E
k_release = 1.185
# Rm-cPAcP	→	Ar-cPAcP
k_elongate = 1.31
# Ar-cPAcP	→	cPAcP+40S+60S+mRNA-cPAcP
k_terminate = 1.612
# cPAcP	→	[]
k_prot_deg = 1.208e-2
create_translation_rules(cPAcP, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)

# *** sPAcP transcription ***
# g-sPAcP+RNAp	↔	g-sPAcP-RNAp
# g-sPAcP-RNAp	→	g-sPAcP+RNAp+mRNA-sPAcP
kf_kr_kcat = [[4.894e-2, 1.238e-3, 0.6989]]  # basal
# mRNA-sPAcP	→	[]
k_mRNA_deg = 0.112
# AR-p-2+g-sPAcP	↔	AR-p-2-g-sPAcP
# AR-p-DHT-2+g-sPAcP	↔	AR-p-DHT-2-g-sPAcP
# AR-p-DHT-AR-p+g-sPAcP	↔	AR-p-DHT-AR-p-g-sPAcP
# AR-p-T-2+g-sPAcP	↔	AR-p-T-2-g-sPAcP
# AR-p-AR-p-T+g-sPAcP	↔	AR-p-AR-p-T-g-sPAcP
# AR-p-DHT-AR-p-T+g-sPAcP	↔	AR-p-DHT-AR-p-T-g-sPAcP
k_tf_on_off = [[9.585e-2, 2.558e-2], [9.074e-6, 0.1017], [5.658e-6, 0.1447], [4.345e-3, 7.3e-2],
               [1.811e-2, 6.045e-3], [1.125e-2, 2.44e-2]]
create_transcription_rules(sPAcP, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species_AR, tf_names_AR)),
                           k_tf_on_off=k_tf_on_off)

# *** sPAcP translation ***
# mRNA-sPAcP+eIF4E	↔	mRNA-sPAcP-eIF4E
# mRNA-sPAcP-eIF4E+40S	↔	mRNA-sPAcP-eIF4E-40S
# mRNA-sPAcP-eIF4E-40S+60S	↔	mRNA-sPAcP-eIF4E-40S-60S
kf_kr = [
    [3.473e-4, 2.741e-2],  # eIF4E binding
    [0.3844, 7.084e-3],     # 40S binding
    [1.532e-2, 2.415e-3]     # 60S binding
]
# mRNA-sPAcP-eIF4E-40S-60S	→	Rm-sPAcP+eIF4E
k_release = 1.089
# Rm-sPAcP	→	Ar-sPAcP
k_elongate = 1.572
# Ar-sPAcP	→	sPAcP+40S+60S+mRNA-sPAcP
k_terminate = 1.474
# sPAcP	→	[]
k_prot_deg = 7.162e-4
create_translation_rules(sPAcP, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)

# *** CycD transcription ***
# g-CycD+RNAp	↔	g-CycD-RNAp
# g-CycD-RNAp	→	g-CycD+RNAp+mRNA-CycD
kf_kr_kcat = [[4.952e-5, 0.1104, 1.099e-2]]  # basal
# mRNA-CycD	→	[]
k_mRNA_deg = 0.8094
# g-CycD+ETS-p	↔	g-CycD-ETS-p
# g-CycD+AP1-p	↔	g-CycD-AP1-p
tf_names_CycD = ['ETS_p', 'AP1_p']
tf_species_CycD = [ETS(erk=None,gene=None,state='p'), AP1(erk=None,gene=None,state='p')]
k_tf_on_off = [[0.138, 1.26], [0.3726, 2.171]]  # ETS-p, AP1-p
# g-CycD-ETS-p+RNAp	↔	g-CycD-ETS-p-RNAp
# g-CycD-ETS-p-RNAp	→	g-CycD-ETS-p+RNAp+mRNA-CycD
# g-CycD-AP1-p+RNAp	↔	g-CycD-AP1-p-RNAp
# g-CycD-AP1-p-RNAp	→	g-CycD-AP1-p+RNAp+mRNA-CycD
kf_kr_kcat.extend(
    [[0.2931, 9.037e-3, 1.156e-2],  # ETS-p
     [0.4288, 2.945e-2, 3.292e-2]]  # AP1-p
)
create_transcription_rules(CycD, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species_CycD, tf_names_CycD)),
                           k_tf_on_off=k_tf_on_off)

# *** CycD translation ***
# mRNA-CycD+eIF4E	↔	mRNA-CycD-eIF4E
# mRNA-CycD-eIF4E+40S	↔	mRNA-CycD-eIF4E-40S
# mRNA-CycD-eIF4E-40S+60S	↔	mRNA-CycD-eIF4E-40S-60S
kf_kr = [
    [2.137e-2, 8.14e-3],   # eIF4E
    [8.773e-2, 5.877e-3],  # 40S
    [0.7418, 1.306e-3]     # 60S
]
# mRNA-CycD-eIF4E-40S-60S	→	Rm-CycD+eIF4E
k_release = 1.194
# Rm-CycD	→	Ar-CycD
k_elongate = 0.9399
# Ar-CycD	→	CycD+40S+60S+mRNA-CycD
k_terminate = 1.781
# CycD	→	[]
k_prot_deg = 6.123e-3
create_translation_rules(CycD, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)

# *** PSA transcription ***
# g-PSA+RNAp	↔	g-PSA-RNAp
# g-PSA-RNAp	→	g-PSA+RNAp+mRNA-PSA
kf_kr_kcat = [[9.158e-9, 1.29e-4, 2.873e-4]]  # basal
# mRNA-PSA	→	[]
k_mRNA_deg = 0.1389
# g-PSA+AR-p-2	↔	g-PSA-AR-p-2
# g-PSA+AR-p-DHT-2	↔	g-PSA-AR-p-DHT-2
# g-PSA+AR-p-DHT-AR-p	↔	g-PSA-AR-p-DHT-AR-p
# g-PSA+AR-p-T-2	↔	g-PSA-AR-p-T-2
# g-PSA+AR-p-AR-p-T	↔	g-PSA-AR-p-AR-p-T
# g-PSA+AR-p-DHT-AR-p-T	↔	g-PSA-AR-p-DHT-AR-p-T
k_tf_on_off = [[0.1372, 4.913e-3], [8.626e-2, 2.007e-4], [1.025e-2, 1.353e-4], [9.22e-5, 8.384e-4],
               [4.169e-5, 5.429e-3], [2.392e-4, 5.202e-4]]
# g-PSA-AR-p-2+RNAp	↔	g-PSA-AR-p-2-RNAp
# g-PSA-AR-p-2-RNAp	→	g-PSA-AR-p-2+RNAp+mRNA-PSA
# g-PSA-AR-p-DHT-2+RNAp	↔	g-PSA-AR-p-DHT-2-RNAp
# g-PSA-AR-p-DHT-2-RNAp	→	g-PSA-AR-p-DHT-2+RNAp+mRNA-PSA
# g-PSA-AR-p-DHT-AR-p+RNAp	↔	g-PSA-AR-p-DHT-AR-p-RNAp
# g-PSA-AR-p-DHT-AR-p-RNAp	→	g-PSA-AR-p-DHT-AR-p+RNAp+mRNA-PSA
# g-PSA-AR-p-T-2+RNAp	↔	g-PSA-AR-p-T-2-RNAp
# g-PSA-AR-p-T-2-RNAp	→	g-PSA-AR-p-T-2+RNAp+mRNA-PSA
# g-PSA-AR-p-AR-p-T+RNAp	↔	g-PSA-AR-p-AR-p-T-RNAp
# g-PSA-AR-p-AR-p-T-RNAp	→	g-PSA-AR-p-AR-p-T+RNAp+mRNA-PSA
# g-PSA-AR-p-DHT-AR-p-T+RNAp	↔	g-PSA-AR-p-DHT-AR-p-T-RNAp
# g-PSA-AR-p-DHT-AR-p-T-RNAp	→	g-PSA-AR-p-DHT-AR-p-T+RNAp+mRNA-PSA
kf_kr_kcat.extend(
    [[1.182e-4, 8.565e-5, 8.167e-2],  # AR-p-2
     [7.817e-2, 1.577e-4, 2.238e-2],  # AR-p-DHT-2
     [1.091e-4, 1.029e-5, 3.064e-4],  # AR-p-DHT-AR-p
     [4.362e-4, 2.375e-4, 8.891e-5],  # AR-p-T-2
     [8.539e-5, 1.204e-4, 1.959e-4],  # AR-p-AR-p-T
     [2.092e-4, 1.516e-4, 2.128e-4]]  # R-p-DHT-AR-p-T
)
create_transcription_rules(PSA, kf_kr_kcat, k_mRNA_deg, tfs=list(zip(tf_species_AR, tf_names_AR)),
                           k_tf_on_off=k_tf_on_off)

# *** PSA translation ***
# mRNA-PSA+eIF4E	↔	mRNA-PSA-eIF4E
# mRNA-PSA-eIF4E+40S	↔	mRNA-PSA-eIF4E-40S
# mRNA-PSA-eIF4E-40S+60S	↔	mRNA-PSA-eIF4E-40S-60S
kf_kr = [
    [1.278e-4, 5.894e-6],  # eIF4E
    [9.861e-4, 4.486e-6],  # 40S
    [7.736e-4, 6.108e-4]   # 60S
]
# mRNA-PSA-eIF4E-40S-60S	→	Rm-PSA+eIF4E
k_release = 3.98e-4
# Rm-PSA	→	Ar-PSA
k_elongate = 7.268e-3
# Ar-PSA	→	PSA+40S+60S+mRNA-PSA
k_terminate = 3.797e-2
# PSA	→	[]
k_prot_deg = 9.309e-6
create_translation_rules(PSA, kf_kr, k_release, k_elongate, k_terminate, k_prot_deg)


# === OBSERVABLES ===

Observable('Lig_free', EGF(r=None))
Observable('Lig_bound', EGF(r=ANY))
Observable('Rec_unphos', EGFR(state='u'))
Observable('Rec_phos', EGFR(state='p'))

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from pysb.simulator import ScipyOdeSimulator

    # simulation commands
    tspan = np.linspace(0, 10, 101)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    result = sim.run()

    rules_unidirectional = [rule for rule in model.rules if rule.rate_reverse is None]
    rules_bidirectional = [rule for rule in model.rules if rule.rate_reverse is not None]
    print('# of rules:', len(rules_unidirectional) + 2 * len(rules_bidirectional))
    print('# of reactions:', len(model.reactions))

    rule_names = np.array([rule.name for rule in model.rules])
    rxn_rules = np.unique([rule_name for rxn in model.reactions for rule_name in rxn['rule']])
    print(len(set(rule_names)-set(rxn_rules)))
    print(np.array(list(set(rule_names) - set(rxn_rules))))


    # plot results
    plt.figure(constrained_layout=True)
    for obs in model.observables:
        plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.legend(loc='best')

    plt.show()
