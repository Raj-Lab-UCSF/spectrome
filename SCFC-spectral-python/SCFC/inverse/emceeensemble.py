import emcee #the affine invariant MCMC sampler from Foreman-Mackey, D., Hogg, D. W., Lang, D., Goodman, J., Feb. 2012. emcee: The MCMC Hammer.
import bayesian_funcs as bayes

#this is far from functioning yet, I am just trying to get all the pieces together
#to start making an outline of how it will all fit, and build it up from there.

class EmceeWalkers():

    def __init__(self,ndim, nwalkers, nsteps):
        self.ndim = ndim
        self.nwalkers = nwalkers
        self.nsteps = nsteps

    def set_fixed_params(self):
        #here we are going to set up stuff that will be fixed in our model as we walk through
        #parameter spaceself.

    def start_position(self):
        #find start point of nwalkers

    def set_ensemble(self):
        #set up the emcee ensemble sampler by passing the functions, additional stuff
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, bayes.lnprob, args=(x, y, yerr))

    def run_ensemble(self):
        pos = self.start_position()
        sampler.run_mcmc(pos, self.nsteps)
