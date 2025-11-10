from pysb import *
from pysb import MonomerPattern
from pysb.util import alias_model_components
import numpy as np

model = None


def set_model(m):
    global model
    model = m


def create_transcription_rules(prot_monomer, kf_kr_kcat, k_deg, tfs=None, k_tf_on_off=None):
    alias_model_components()

    # Transcription Scheme:
    # reaction			                                                description
    # g-ProteinX + RNAp	    <->	g-ProteinX-RNAp	                        RNA polymerase binds to ProteinX gene
    # g-ProteinX-RNAp	    ->	g-ProteinX+RNAp + mRNA-ProteinX	        transcription of ProteinX mRNA
    # g-ProteinX + TF       <->	g-ProteinX-TF	                        transcription factor (TF) binds to ProteinX gene
    # g-ProteinX-TF + RNAp	<->	g-ProteinX-TF-RNAp	                    TF-regulated binding of RNA polymerase
    # g-ProteinX-TF-RNAp	->	g-ProteinX-TF + RNAp + mRNA-ProteinX	formation of ProteinX mRNA
    # mRNA-ProteinX         ->  None                                    ProteinX mRNA degradation

    if len(np.array(kf_kr_kcat).shape) == 1:
        kf_kr_kcat = [kf_kr_kcat]

    if len(np.array(tfs).shape) == 0:
        tfs = [tfs]
        k_tf_on_off = [k_tf_on_off]

    gene = model.monomers['g_%s' % prot_monomer.name]
    mrna = model.monomers['mRNA_%s' % prot_monomer.name]

    # RNAp binds to gene and transcribes mRNA
    kf, kr, kcat = \
        [Parameter('%s_%s_RNAp' % (k, gene.name), kf_kr_kcat[0][i]) for i, k in enumerate(['kf', 'kr', 'kcat'])]
    Rule('%s_binds_RNAp' % gene.name,
         gene(rnap=None, tf=None) + RNAp(gene=None) |
         gene(rnap=1, tf=None) % RNAp(gene=1), kf, kr)
    Rule('%s_RNAp_transcribes' % gene.name,
         gene(rnap=1, tf=None) % RNAp(gene=1) >>
         gene(rnap=1, tf=None) % RNAp(gene=1) + mrna(eif4e=None,_40s=None,elong='i'), kcat)

    # mRNA degradation
    k_deg = Parameter('kdeg_%s' % mrna.name, k_deg)
    Rule('%s_degrades' % mrna.name, mrna(eif4e=None,_40s=None,elong='i') >> None, k_deg)

    if tfs is not None:
        for i, tf in enumerate(tfs):
            # transcription factor (TF) binds to gene
            k_on, k_off = [Parameter('%s_%s_%s' % (k, gene.name, tf.monomer.name), k_tf_on_off[i][j])
                           for j, k in enumerate(['kon', 'koff'])]
            Rule('%s_binds_%s' % (gene.name, tf.monomer.name),
                 gene(rnap=None, tf=None) + tf(gene=None) |
                 gene(rnap=None, tf=1) % tf(gene=1), k_on, k_off)

            # RNAp binds to TF-bound gene and transcribes mRNA
            kf, kr, kcat = [Parameter('%s_%s_%s_RNAp' % (k, gene.name, tf.monomer.name), kf_kr_kcat[i+1][j])
                                      for j, k in enumerate(['kf', 'kr', 'kcat'])]
            Rule('%s_%s_binds_RNAp' % (gene.name, tf.monomer.name),
                 gene(rnap=None, tf=1) % tf(gene=1) + RNAp(gene=None) |
                 gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2), kf, kr)
            Rule('%s_%s_RNAp_transcribes' % (gene.name, tf.monomer.name),
                 gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2) >>
                 gene(rnap=2, tf=1) % tf(gene=1) % RNAp(gene=2) + mrna(eif4e=None, _40s=None,elong='i'), kcat)


def create_translation_rules(prot_monomer, kf_kr, k_release, k_elongate, k_terminate, k_deg):
    alias_model_components()

    # Translation Scheme:
    # reaction			                                            description
    # mRNA-ProteinX + eIF4E	        <->	mRNA-ProteinX-eIF4E	        ProteinX mRNA binds initiation factor eIF4E
    # mRNA-ProteinX-eIF4E + 40S	    <->	mRNA-ProteinX-eIF4E-40S	    ProteinX mRNA binds ribosome 40S
    # mRNA-ProteinX-eIF4E-40S + 60S	<->	mRNA-ProteinX-eIF4E-40S-60S	ProteinX mRNA binds ribosome 60S
    # mRNA-ProteinX-eIF4E-40S-60S	->	Rm-ProteinX + eIF4E	        ProteinX mRNA releases initiation factor eIF4E
    # Rm-ProteinX   ->	Ar-ProteinX	                                Elongation phase of translation
    # Ar-ProteinX   ->	ProteinX + 40S + 60S + mRNA-ProteinX	    Termination
    # ProteinX      ->  None                                        Protein degradation

    # Monomer('eIF4E', ['mrna_4ebp1'])
    # Monomer('_40S', ['mrna', '_60s'])
    # Monomer('_60S', ['_40s'])

    mrna = model.monomers['mRNA_%s' % prot_monomer.name]

    site_conditions = dict([(site, None) if site not in prot_monomer.site_states.keys() else
                            (site, prot_monomer.site_states[site][0]) for site in prot_monomer.sites])
    protein = MonomerPattern(prot_monomer, site_conditions, None)

    kf, kr = [Parameter('%s_%s_binds_eIF4E' % (k, mrna.name), val) for k, val in zip(['kf', 'kr'], kf_kr[0])]
    Rule('%s_binds_eIF4E' % mrna.name,
         mrna(eif4e=None,_40s=None,elong='i') + eIF4E(mrna_4ebp1=None) |
         mrna(eif4e=1,_40s=None,elong='i') % eIF4E(mrna_4ebp1=1), kf, kr)

    kf, kr = [Parameter('%s_%s_eIF4E_binds_40S' % (k, mrna.name), val) for k, val in zip(['kf', 'kr'], kf_kr[1])]
    Rule('%s_eIF4E_binds_40S' % mrna.name,
         mrna(eif4e=1,_40s=None, elong='i') % eIF4E(mrna_4ebp1=1) + _40S(mrna=None,_60s=None) |
         mrna(eif4e=1,_40s=2, elong='i') % eIF4E(mrna_4ebp1=1) % _40S(mrna=2,_60s=None), kf, kr)

    kf, kr = [Parameter('%s_%s_eIF4E_40S_binds_60S' % (k, mrna.name), val) for k, val in zip(['kf', 'kr'], kf_kr[2])]
    Rule('%s_eIF4E_40S_binds_60S' % mrna.name,
         mrna(eif4e=1,_40s=2, elong='i') % eIF4E(mrna_4ebp1=1) % _40S(mrna=2, _60s=None) + _60S(_40s=None) |
         mrna(eif4e=1,_40s=2, elong='i') % eIF4E(mrna_4ebp1=1) % _40S(mrna=2, _60s=3) % _60S(_40s=3), kf, kr)

    k_release = Parameter('k_release_%s' % prot_monomer.name, k_release)
    Rule('%s_eIF4E_40S_60S_releases_eIF4E' % mrna.name,
         mrna(eif4e=1,_40s=2, elong='i') % _40S(mrna=2, _60s=3) % _60S(_40s=3) % eIF4E(mrna_4ebp1=1) >>
         mrna(eif4e=None,_40s=2, elong='i') % _40S(mrna=2, _60s=3) % _60S(_40s=3) + eIF4E(mrna_4ebp1=None) , k_release)

    k_elongate = Parameter('k_elongate_%s' % prot_monomer.name, k_elongate)
    Rule('%s_40S_60S_elongates' % mrna.name,
         mrna(eif4e=None, _40s=2, elong='i') % _40S(mrna=2, _60s=3) % _60S(_40s=3) >>
         mrna(eif4e=None, _40s=2, elong='a') % _40S(mrna=2, _60s=3) % _60S(_40s=3), k_elongate)

    k_terminate = Parameter('k_terminate_%s' % prot_monomer.name, k_terminate)
    Rule('%s_40S_60S_translates' % mrna.name,
         mrna(eif4e=None, _40s=2, elong='a') % _40S(mrna=2, _60s=3) % _60S(_40s=3) >>
         mrna(eif4e=None, _40s=None, elong='i') + _40S(mrna=None, _60s=None) + _60S(_40s=None) + protein, k_terminate)

    k_deg = Parameter('k_deg_%s' % prot_monomer.name, k_deg)
    Rule('%s_degrades' % prot_monomer.name, protein >> None, k_deg)