from pysb import *
from pysb import MonomerPattern, ComplexPattern
from pysb.util import alias_model_components
import numpy as np
import re

model = None


def set_model(m):
    global model
    model = m


def _get_gene_unbound_bound_tf(gene, tf):
    if isinstance(tf, MonomerPattern):
        gene_unbound_tf = gene(rnap=None, tf=None)
        gene_bound_tf = gene(rnap=None, tf=50)
        tf_bound = tf(gene=50)
    elif isinstance(tf, ComplexPattern):
        site_conditions = dict(
            [('rnap', None)] + [('tf%d' % (i + 1), None) for i in range(len(gene.sites) - 1)])
        gene_unbound_tf = MonomerPattern(gene, site_conditions, None)
        site_conditions = dict(
            [('rnap', None)] + [('tf%d' % (i + 1), 50 + i) for i in range(len(gene.sites) - 1)])
        gene_bound_tf = MonomerPattern(gene, site_conditions, None)
        mp_gene_sites = [mp for mp in tf.monomer_patterns if 'gene' in mp.monomer.sites]
        tf_bound = ComplexPattern([mp(gene=50 + i) for i, mp in enumerate(mp_gene_sites)] +
                                  [mp for mp in tf.monomer_patterns if mp not in mp_gene_sites],
                                  None)
    return gene_unbound_tf, gene_bound_tf, tf_bound


def create_transcription_rules(prot_monomer, kf_kr_kcat, k_deg, tfs=None, k_tf_on_off=None):
    alias_model_components()

    # Transcription Scheme:
    # reaction			                                                description
    # g-ProteinX + RNAp	    <->	g-ProteinX-RNAp	                        RNA polymerase binds to ProteinX gene
    # g-ProteinX-RNAp	    ->	g-ProteinX + RNAp + mRNA-ProteinX	    transcription of ProteinX mRNA
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
    tf_sites = [site for site in gene.sites if 'tf' in site]
    tf_sites_dict = dict(list(zip(tf_sites, [None] * len(tf_sites))))
    # g-ProteinX + RNAp	    <->	g-ProteinX-RNAp	                        RNA polymerase binds to ProteinX gene
    Rule('%s_binds_RNAp' % gene.name,
         gene(rnap=None, **tf_sites_dict) + RNAp(gene=None) |
         gene(rnap=1, **tf_sites_dict) % RNAp(gene=1), kf, kr)
    # g-ProteinX-RNAp	    ->	g-ProteinX + RNAp + mRNA-ProteinX	    transcription of ProteinX mRNA
    Rule('%s_RNAp_transcribes' % gene.name,
         gene(rnap=1, **tf_sites_dict) % RNAp(gene=1) >>
         gene(rnap=None, **tf_sites_dict) + RNAp(gene=None) + mrna(eif4e=None,_40s=None,elong='i'), kcat)

    # mRNA degradation
    k_deg = Parameter('kdeg_%s' % mrna.name, k_deg)
    # mRNA-ProteinX         ->  None                                    ProteinX mRNA degradation
    Rule('%s_degrades' % mrna.name, mrna(eif4e=None,_40s=None,elong='i') >> None, k_deg)

    if tfs is not None:
        for i, (tf, tf_name) in enumerate(tfs):
            # transcription factor (TF) binds to gene
            k_on, k_off = [Parameter('%s_%s_%s' % (k, gene.name, tf_name), k_tf_on_off[i][j])
                           for j, k in enumerate(['kon', 'koff'])]
            gene_unbound_tf, gene_bound_tf, tf_bound = _get_gene_unbound_bound_tf(gene, tf)
            # g-ProteinX + TF       <->	g-ProteinX-TF	                        transcription factor (TF) binds to ProteinX gene
            Rule('%s_binds_%s' % (gene.name, tf_name),
                 gene_unbound_tf + tf | gene_bound_tf % tf_bound, k_on, k_off)

            # RNAp binds to TF-bound gene and transcribes mRNA
            if len(kf_kr_kcat) > i + 1:
                kf, kr, kcat = [Parameter('%s_%s_%s_RNAp' % (k, gene.name, tf_name), kf_kr_kcat[i+1][j])
                                          for j, k in enumerate(['kf', 'kr', 'kcat'])]
                # g-ProteinX-TF + RNAp	<->	g-ProteinX-TF-RNAp	                    TF-regulated binding of RNA polymerase
                Rule('%s_%s_binds_RNAp' % (gene.name, tf_name),
                     gene_bound_tf % tf_bound + RNAp(gene=None) |
                     gene_bound_tf(rnap=100) % tf_bound % RNAp(gene=100), kf, kr)
                # g-ProteinX-TF-RNAp	->	g-ProteinX-TF + RNAp + mRNA-ProteinX	formation of ProteinX mRNA
                Rule('%s_%s_RNAp_transcribes' % (gene.name, tf_name),
                     gene_bound_tf(rnap=100) % tf_bound % RNAp(gene=100) >>
                     gene_bound_tf(rnap=None) % tf_bound + RNAp(gene=None) + mrna(eif4e=None, _40s=None, elong='i'),
                     kcat)


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


def _extract_bng_multipliers(mult_factor_lines):
    factor_dict = {}
    for line in mult_factor_lines.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        rate_expr = fields[3]  # e.g. "0.5*kf_EGFR_EGF_dimerize"
        m = re.match(r"^([0-9.]+)\*([A-Za-z_]\w*)$", rate_expr)
        if not m:
            continue
        factor = float(m.group(1))
        par_name = m.group(2)
        factor_dict[par_name] = factor
    return factor_dict


def divide_out_bng_multipliers(model, mult_factor_lines, verbose=True):
    factor_dict = _extract_bng_multipliers(mult_factor_lines)
    par_names = [p.name for p in model.parameters]
    for par_name, factor in factor_dict.items():
        if par_name not in par_names:
            raise Exception(f"Parameter not found: {par_name}")
        old_value = model.parameters[par_name].value
        new_value = old_value / factor
        model.parameters[par_name].value = new_value
        if verbose:
            print(f"{par_name}: {old_value} / {factor} = {new_value}")
    return factor_dict


def _parse_flat_reactions_from_source(filename):


    def _split_flat_side(side):
        """
        Convert:
            '2*EGFR-EGF' -> ['EGFR-EGF', 'EGFR-EGF']
            '4*Her2-2-p+cPAcP-4' -> ['Her2-2-p', 'Her2-2-p', 'Her2-2-p', 'Her2-2-p', 'cPAcP-4']
        """
        side = side.strip()
        if side == "[]":
            return []
        species = []
        for token in side.split("+"):
            token = token.strip()
            if not token or token == "[]":
                continue
            m = re.match(r"^(\d+)\*(.+)$", token)
            if m:
                coeff = int(m.group(1))
                name = m.group(2).strip()
                species.extend([name] * coeff)
            else:
                species.append(token)
        return species


    flat_reactions = []
    arrow_pattern = re.compile(r"#\s*(\d+)\s*[.\t ]+(.*?)\s*(↔|->|→)\s*(.*)")
    with open(filename, "r") as f:
        for line in f:
            m = arrow_pattern.match(line)
            if not m:
                continue
            idx, lhs, arrow, rhs = m.groups()
            idx = int(idx)  # convert to integer
            # remove trailing kinetic info after the actual reaction
            rhs = re.split(r"\s{2,}|\t", rhs.strip())[0]
            flat_reactions.append({
                "index": idx,
                "reactants": _split_flat_side(lhs),
                "products": _split_flat_side(rhs),
                "arrow": arrow,
            })
    return flat_reactions


def _pattern_in_species_dict(pattern, species_dict):
    for key in species_dict.keys():
        if pattern.is_equivalent_to(key):
            return key
    return False


if __name__ == "__main__":
    from AR_model import model

    flat_reactions = _parse_flat_reactions_from_source('AR_model.py')

    species_dict = {}
    conflicts = []

    for flat_rxn, rule in zip(flat_reactions, model.rules):
        print(flat_rxn)
        print(rule)
        flat_species = flat_rxn["reactants"] + flat_rxn["products"]

        pysb_patterns = (
                list(rule.reactant_pattern.complex_patterns) +
                list(rule.product_pattern.complex_patterns)
        )

        if len(flat_species) != len(pysb_patterns):
            msg = f"Length mismatch for rule {rule.name}\n" + "Flat: " + str(flat_species) + \
                  "\nPySB:" + str(pysb_patterns)
            raise Exception(msg)

        for pattern, name in zip(pysb_patterns, flat_species):
            print('%s: %s' % (pattern, name))
            key = _pattern_in_species_dict(pattern, species_dict)
            print('key:', key)
            if key is False:
                species_dict[pattern] = name
                print('ADDED %s' % name)
            elif species_dict[key] != name:
                # conflicts.append((key, species_dict[key], name))
                raise Exception('Pattern %s associated with multiple names (%s, %s)' %
                                (pattern, species_dict[key], name))
        print()

    # print the dictionary assignments
    for i, (pattern, name) in enumerate(species_dict.items()):
        print(f"{i}: species_dict[{pattern}] = {name!r}")

    # quit()

    #####
    from pysb import as_complex_pattern
    from pysb.simulator import ScipyOdeSimulator
    sim = ScipyOdeSimulator(model, verbose=True, cleanup=True)
    remove_idxs = []
    for i, species in enumerate(model.species):
        found = False
        for pattern in species_dict.keys():
            if species.is_equivalent_to(pattern):
                found = True
                break
        if not found:
            print('Removing:')
            print('   %d: %s (%s = %g)' % \
                  (i, model.initials[i].pattern, model.initials[i].value.name, model.initials[i].value.value))
            remove_idxs.append(i)
    model.species = [v for i, v in enumerate(model.species) if i not in remove_idxs]

    print()
    for i, species in enumerate(model.species):
        for key in species_dict.keys():
            if species.is_equivalent_to(key):
                print('%d:' % i, species_dict[key])
                break
    print()
    print('Number of species:', len(model.species))
    print()
    # for sp in model.species:
    #     print(sp)

    quit()

    i = 0
    for rule in model.rules:
        idx = str(i)
        i += 1
        if rule.is_reversible:
            idx += ',%d' % i
            i += 1
        # print('%s:' % idx, rule)
        # print('Reactants:')
        reactants = []
        for cp in rule.reactant_pattern.complex_patterns:
            # print(cp, type(cp))
            # print(cp.is_concrete())
            for key in species_dict.keys():
                if as_complex_pattern(key).is_equivalent_to(cp):
                    reactants.append(species_dict[key])
                    break
            # print(reactants)
        # print('Products:')
        products = []
        for cp in rule.product_pattern.complex_patterns:
            # print(cp, type(cp))
            # print(cp.is_concrete())
            for key in species_dict.keys():
                if as_complex_pattern(key).is_equivalent_to(cp):
                    products.append(species_dict[key])
                    break
            # print(products)
        # print rule to screen
        arrow = '<->' if rule.is_reversible else '-->'
        rxn_string = ' + '.join(reactants) + ' %s ' % arrow + ' + '.join(products)
        print('%s:' % idx, rxn_string)
        # print()
