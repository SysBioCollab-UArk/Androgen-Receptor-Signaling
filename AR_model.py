from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

# monomers
Monomer('EGF',['r'])
Monomer('EGFR',['l','d','state'],{'state':['u','p']})


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
