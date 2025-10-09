from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

# monomers

# Ligands & Receptors
Monomer('EGF',['r'])
Monomer('EGFR',['l','d','state'],{'state':['u','p']})

# Androgen axis (Ligand + receptor + co-regulators)
Monomer('DHT', ['r'])
Monomer('AR', ['l'])
Monomer('PGC1A', ['b'])
Monomer('ESRRA', ['b'])
Monomer('NRF1', ['b'])
Monomer('TFAM', ['b'])

# PI3K-AKT axis and crosstalk
Monomer('PI3K', ['b'])
Monomer('PIP3', ['b'])
Monomer('AKT', ['b'])
Monomer('PTEN', ['b'])
Monomer('mTOR', ['b'])
Monomer('FOXO', ['b'])
Monomer('GSK3B', ['b'])
Monomer('BAD', ['b'])

# MAPK axis (optional downstream of EGFR/RAS)
Monomer('GRB2', ['b'])
Monomer('SOS', ['b'])
Monomer('RAS', ['b'])
Monomer('RAF', ['b'])
Monomer('MEK', ['b'])
Monomer('ERK', ['b'])

# Mitochondrial dynamics & quality control
Monomer('DRP1', ['b'])
Monomer('MFN1', ['b'])
Monomer('MFN2', ['b'])
Monomer('OPA1', ['b'])
Monomer('PINK1', ['b'])
Monomer('PRKN', ['b'])
Monomer('ROS', ['b'])

# Mitochondrial UPR (UPRmt)
Monomer('HSPD1', ['b'])
Monomer('HSPE1', ['b'])
Monomer('HSPA9', ['b'])
Monomer('LONP1', ['b'])
Monomer('CLPP', ['b'])
Monomer('ATF5', ['b'])

# Apoptosis (intrinsic pathway)
Monomer('CYCS', ['b'])
Monomer('APAF1', ['b'])
Monomer('CASP9', ['b'])
Monomer('CASP3', ['b'])
Monomer('XIAP', ['b'])
Monomer('SMAC', ['b'])

# BCL-2 family (anti-apoptotic, effectors, BH3-only)
Monomer('BCL2', ['bh3'])
Monomer('BCLXL', ['bh3'])
Monomer('MCL1', ['bh3'])
Monomer('BAX', ['b'])
Monomer('BAK', ['b'])
Monomer('BIM', ['bh3'])
Monomer('PUMA', ['bh3'])
Monomer('NOXA', ['bh3'])


# initial concentrations
Parameter('Lig_0', 100)
Parameter('Rec_0', 100)

Initial(EGF(r=None), Lig_0)
Initial(EGFR(l=None, d=None, state='u'), Rec_0)

# rules
Parameter('kf_LR', 1)
Parameter('kr_LR', 1)
Rule('lig_binds_rec', EGF(r=None) + EGFR(l=None, d=None) | EGF(r=1) % EGFR(l=1, d=None), kf_LR, kr_LR)


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
