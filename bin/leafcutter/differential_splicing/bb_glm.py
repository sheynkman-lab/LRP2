import torch
import torch.nn.functional as F

import pyro
import pyro.distributions as dist
from pyro.infer.autoguide import AutoDiagonalNormal, AutoGuideList, AutoDelta, init_to_value
from pyro.infer import SVI, Trace_ELBO, config_enumerate, infer_discrete
from pyro.ops.indexing import Vindex

from pyro import poutine

from torch.distributions import constraints

import numpy as np
import scipy.stats
import time
from sklearn import linear_model

from leafcutter.differential_splicing.optim import fit_with_lbfgs, fit_with_SVI

from dataclasses import dataclass

#torch.set_num_threads(20)

@dataclass
class LeafcutterFit:
    """
    Stores the result of fitting the LeafCutter differential splicing model to all junctions. 
    """
    beta: torch.Tensor # length J vector
    conc: torch.Tensor # length J vector
    loss: float # will want this per junction somehow
    likelihoods: torch.Tensor 
    exit_status: str

def convertr(hyperparam, name, device): 
    return pyro.sample(name, hyperparam) if (type(hyperparam) != float) else torch.tensor(hyperparam, device = device) 

class BayesianBetaBinomialModel(pyro.nn.PyroModule):
    """
    The LeafCutter differential splicing model for one cluster. 
    """
    
    def __init__(self, P_null, P_full, J, eps = 1e-8, gamma_shape = 1.0001, gamma_rate = 1e-4, beta_scale = np.inf, prior_prob = 0.1, multiconc = True): 
        """
        Initialize the BetaBinomialModel.

        Args:
            P (int): Number of covariates.
            J (int): Number of junctions.
            eps (float): A small constant to prevent numerical issues.
            gamma_shape (float): Shape parameter for the gamma distribution used for the concentration prior. 
            gamma_rate (float): Rate parameter for the gamma distribution used for the concentration prior. 
            beta_scale (float): prior std for coefficients. 
            multiconc (bool): Indicates whether to use a separate concentration parameter for each junction.

        """
        super().__init__()
        self.P_full = P_full
        self.P_null = P_null
        self.J = J
        self.eps = eps
        self.binomial = gamma_shape is None
        self.gamma_shape = gamma_shape
        self.gamma_rate = gamma_rate
        self.multiconc = multiconc
        self.beta_full_scale = beta_scale
        self.beta_null_scale = beta_scale
        self.prior_prob = prior_prob

    @config_enumerate
    def forward(self, x_null, x_full, y, n): 
        """
        Run the generative process. 

        Args:
            x_null (torch.Tensor): Input data representing covariates. N x P_null.
            x_full (torch.Tensor): Input data representing covariates. N x P_full.
            y (torch.Tensor): Observed data, i.e. junction counts. N x J.  
            n (torch.Tensor): ATSE (total) counts. N x J. 
        """

        device = x_full.device

        beta_null_scale = convertr(self.beta_null_scale, "beta_null_scale", device = device)
        with pyro.plate("P_null", self.P_null): # beta is P (covariates) x J (junctions)
            beta_dist = dist.Normal(0.0, beta_null_scale) if torch.isfinite(beta_null_scale) else dist.ImproperUniform(constraints.real, (), ())
            beta_null = pyro.sample("beta_null", beta_dist.expand([self.P_null, self.J]).to_event(1)) 

        beta_full_scale = convertr(self.beta_full_scale, "beta_full_scale", device = device)
        with pyro.plate("P_full", self.P_full): # beta is P (covariates) x J (junctions)
            beta_dist = dist.Normal(0.0, beta_full_scale) if torch.isfinite(beta_full_scale) else dist.ImproperUniform(constraints.real, (), ())
            beta_full = pyro.sample("beta_full", beta_dist.expand([self.P_full, self.J]).to_event(1)) 

        logits_full = x_full @ beta_full # logits is N x J
        logits_null = x_null @ beta_null

        prior_prob = convertr(self.prior_prob, "prior_prob", device = device)
        logits_combined = torch.stack((logits_null,logits_full)).transpose(1,2) # 2 x J x N

        # make 2 x J x N, then index to be J x N, then transpose back to N x J
        #logits = logits_combined.transpose(1,2)[assignment.long(), torch.arange(self.J, device=device)].transpose(0,1)

        if self.binomial: 
            likelihood_factors = y * logits - n * F.softplus(logits)
        else:
            gamma_shape = convertr(self.gamma_shape, "gamma_shape", device = device)
            gamma_rate = convertr(self.gamma_rate, "gamma_rate", device = device)
            conc_prior = dist.Gamma(gamma_shape,gamma_rate)
            if self.multiconc: 
                conc_prior = conc_prior.expand([self.J]).to_event(1) 
            conc_param = pyro.sample("conc", conc_prior)
            
            #sum_ab = a + b
            #likelihood_factors = sum_ab.lgamma() + (a+y).lgamma() + (n - y + b).lgamma() - (sum_ab + n).lgamma() - a.lgamma() - b.lgamma()

            
            with pyro.plate("data", self.J): # over junctions
                assignment = pyro.sample(
                    'assignment', 
                    dist.Bernoulli(prior_prob)
                ).long()
                J_arange = torch.arange(self.J, device = device)
                logits = Vindex(logits_combined)[assignment,J_arange,:] # 2 x J x N
                logits = logits.transpose(-1,-2)
                #logits = logits_combined[assignment,J_arange].transpose(0,1)
                p = logits.sigmoid()
                a = p * conc_param + self.eps
                b = (1.-p) * conc_param + self.eps
                y_dist = dist.BetaBinomial(a, b, total_count = n)
                pyro.sample("obs", y_dist, obs=y)

        # manual implementation of the likelihood is faster since it avoids calculating normalization constants
        #pyro.factor("likelihood", likelihood_factors.sum())

        return 1.


class BetaBinomialModel(pyro.nn.PyroModule):
    """
    The LeafCutter differential splicing model for one cluster. 
    """
    
    def __init__(self, P, J, eps = 1e-8, gamma_shape = 1.0001, gamma_rate = 1e-4, beta_scale = np.inf, multiconc = True): 
        """
        Initialize the BetaBinomialModel.

        Args:
            P (int): Number of covariates.
            J (int): Number of junctions.
            eps (float): A small constant to prevent numerical issues.
            gamma_shape (float): Shape parameter for the gamma distribution used for the concentration prior. 
            gamma_rate (float): Rate parameter for the gamma distribution used for the concentration prior. 
            beta_scale (float): prior std for coefficients. 
            multiconc (bool): Indicates whether to use a separate concentration parameter for each junction.

        """
        super().__init__()
        self.P = P
        self.J = J
        self.eps = eps
        self.binomial = gamma_shape is None
        self.gamma_shape = gamma_shape
        self.gamma_rate = gamma_rate
        self.multiconc = multiconc
        self.beta_scale = beta_scale
       
    def forward(self, x, y, n): 
        """
        Run the generative process. 

        Args:
            x (torch.Tensor): Input data representing covariates.
            y (torch.Tensor): Observed data, i.e. junction counts. 
            n (torch.Tensor): ATSE (total) counts. 
        """

        beta_scale = convertr(self.beta_scale, "beta_scale", device = x.device)
        with pyro.plate("covariates", self.P): # beta is P (covariates) x J (junctions)
            beta_dist = dist.Normal(0.0, beta_scale) if torch.isfinite(beta_scale) else dist.ImproperUniform(constraints.real, (), ())
            b_param = pyro.sample("beta", beta_dist.expand([self.P, self.J]).to_event(1)) 

        logits = x @ b_param
        
        if self.binomial: 
            likelihood_factors = y * logits - n * F.softplus(logits)
        else:
            gamma_shape = convertr(self.gamma_shape, "gamma_shape", device = x.device)
            gamma_rate = convertr(self.gamma_rate, "gamma_rate", device = x.device)
            conc_prior = dist.Gamma(gamma_shape,gamma_rate)
            if self.multiconc: 
                conc_prior = conc_prior.expand([self.J]).to_event(1) 
            conc_param = pyro.sample("conc", conc_prior)
            p = logits.sigmoid()
            a = p * conc_param + self.eps
            b = (1.-p) * conc_param + self.eps
            sum_ab = a + b
            
            likelihood_factors = sum_ab.lgamma() + (a+y).lgamma() + (n - y + b).lgamma() - (sum_ab + n).lgamma() - a.lgamma() - b.lgamma()

        # manual implementation of the likelihood is faster since it avoids calculating normalization constants
        pyro.factor("likelihood", likelihood_factors.sum())

        return likelihood_factors.sum(0)

class SimpleGuide():
    """
    Only need a simple guide here, unlike for DM, since there is no weirdness with extra DoF. 
    """
    
    def __init__(self, init_beta, binomial = False, init_conc = None, multiconc = True, conc_max = 300.):
        self.multiconc = multiconc
        self.conc_max = conc_max
        (P,J) = init_beta.shape
        torch_types = { "device" : init_beta.device, "dtype" : init_beta.dtype }
        self.binomial = binomial
        self.init_beta = init_beta
        if not binomial: 
            default_init_conc = torch.full([J],10.,**torch_types) if multiconc else torch.tensor(10.,**torch_types)
            self.init_conc = default_init_conc if (init_conc is None) else init_conc
        
    @property
    def conc(self):
        if self.binomial: return torch.inf
        return pyro.param("conc_loc", lambda: self.init_conc.clone().detach(), constraint=constraints.interval(0., self.conc_max))
    
    @property
    def beta(self):
        return pyro.param("beta_loc", lambda: self.init_beta.clone().detach())
    
    def __call__(self, x, y, n): # ok i think this acutally incorporates the constraint
        if not self.binomial: 
            conc_param = pyro.param("conc_loc", lambda: self.init_conc.clone().detach(), constraint=constraints.interval(0., self.conc_max))
            pyro.sample("conc", dist.Delta(self.conc, event_dim = 1 if self.multiconc else 0)) # do we need to adjust for the constraint manually here? 
        beta_param = self.beta
        with pyro.plate("covariates", beta_param.shape[0]):
            pyro.sample("beta", dist.Delta(beta_param, event_dim = 1))

def brr_initialization(x, y, n, y_eps = 0.33, n_eps = 1.): 
    """
    Try to get a good initialization using Bayesian ridge regression per junction. 
    """
    y_norm = torch.log( (y+y_eps) / (n+n_eps) ).cpu().numpy()
    x_np = x.cpu().numpy()
    
    N,P = x.shape
    J = y_norm.shape[1]
    beta_mm = np.zeros([P,J])
    reg = linear_model.BayesianRidge()
    for j in np.arange(J): 
        reg.fit(x_np, y_norm[:,j])
        beta_mm[:,j] = reg.coef_
    return torch.tensor(beta_mm, dtype = x.dtype, device = x.device)
    
def rr_initialization(x, y, n, regularizer = 0.001, y_eps = 0.33, n_eps = 1.): 
    """
    Try to get a good initialization for beta by moment matching. V slightly slower than BRR overall.
    """
    y_norm = torch.log( (y+y_eps) / (n+n_eps) ).cpu().numpy()
    # get estimate by moment matching
    I = torch.eye(x.shape[1], dtype = x.dtype, device = x.device)
    return torch.linalg.solve( x.t() @ x + regularizer * I, x.t() @ y_norm )
    
def fit_binomial_glm(x, y, n, beta_init = None, fitter = fit_with_lbfgs): 
    """
    Fit a binomial GLM (no overdispersion). 
    """
    [N,P]=x.shape
    J = y.shape[1]
    pyro.clear_param_store()
    
    if beta_init is None: 
        beta_init = torch.zeros(P, J, device = x.device, dtype = x.dtype)

    binomial_model = BetaBinomialModel(P, J, None, gamma_shape = None)
    guide = SimpleGuide(beta_init, binomial = True)
    losses = fitter(binomial_model, guide, [x, y, n])
    return guide.beta

def fit_bb_glm(x, y, n, beta_init, conc_max = 3000., concShape=1.0001, concRate=1e-4, multiconc = True, fitter = fit_with_lbfgs, eps = 1.0e-8): 
    """
    Fits a Beta binomial generalized linear model.

    Args:
        x (torch.Tensor): The input data with shape (N, P) where N is the number of observations and P is the number of predictors.
        y (torch.Tensor): The target data with shape (N, J) where J is the number of splice junctions.
        beta_init (torch.Tensor): The initial values for the coefficients, shape=(P,J)
        conc_max (float, optional): Maximum concentration parameter value. Defaults to 3000.
        concShape (float, optional): Shape parameter for concentration prior. Defaults to 1.0001.
        concRate (float, optional): Rate parameter for concentration prior. Defaults to 1e-4.
        multiconc (bool, optional): Whether to use multiple concentration parameters. Defaults to True.
        fitter (function, optional): The fitting function to use. Defaults to fit_with_lbfgs.
        guide_type (class, optional): The type of guide to use. Defaults to Simpl.
        eps (float): small pseudocount added to DM parameter for numerical stability. 

    Returns:
        LeafcutterFit: An object containing fitted model parameters, losses, and exit status.
    """
    pyro.clear_param_store()

    (N,P) = x.shape
    J = y.shape[1]
    
    model = BetaBinomialModel(P, J, gamma_shape = concShape, gamma_rate = concRate, multiconc = multiconc, eps = eps)
    guide = SimpleGuide(beta_init, multiconc = multiconc, conc_max = conc_max)
    losses, exit_status = fitter(model, guide, [x, y, n])

    # need to get per junction likelihoods
    guide_tr = pyro.poutine.trace(guide).get_trace(x, y, n)
    likelihoods = pyro.poutine.replay(model, trace = guide_tr)(x, y, n)

    return LeafcutterFit(
            beta = guide.beta.clone().detach(), 
            conc = guide.conc.clone().detach(),
            likelihoods = likelihoods,
            loss = losses[-1],
            exit_status = exit_status
        ) 

def svi_bb_glm(x, y, n, beta_init, conc_max = 3000., concShape=1.0001, concRate=1e-4, multiconc = True, beta_scale = np.inf, eps = 1.0e-8): 
    """
    Fits a Beta binomial generalized linear model.

    Args:
        x (torch.Tensor): The input data with shape (N, P) where N is the number of observations and P is the number of predictors.
        y (torch.Tensor): The target data with shape (N, J) where J is the number of splice junctions.
        beta_init (torch.Tensor): The initial values for the coefficients, shape=(P,J)
        conc_max (float, optional): Maximum concentration parameter value. Defaults to 3000.
        concShape (float, optional): Shape parameter for concentration prior. Defaults to 1.0001.
        concRate (float, optional): Rate parameter for concentration prior. Defaults to 1e-4.
        multiconc (bool, optional): Whether to use multiple concentration parameters. Defaults to True.
        fitter (function, optional): The fitting function to use. Defaults to fit_with_lbfgs.
        guide_type (class, optional): The type of guide to use. Defaults to Simpl.
        eps (float): small pseudocount added to DM parameter for numerical stability. 

    Returns:
        LeafcutterFit: An object containing fitted model parameters, losses, and exit status.
    """
    pyro.clear_param_store()

    (N,P) = x.shape
    J = y.shape[1]
    
    model = BetaBinomialModel(P, J, gamma_shape = concShape, gamma_rate = concRate, multiconc = multiconc, beta_scale = beta_scale, eps = eps)
    guide = AutoGuideList(model)
    to_vb = ["conc", "beta"]
    init_dic = { # TODO: GPU 
        "beta": beta_init,
        "conc": torch.full([J], 10.)
    }
    guide.add(AutoDelta(
        poutine.block(model, hide = to_vb),
        init_loc_fn = init_to_value(values=init_dic)))
    guide.add(
        AutoDiagonalNormal(
            poutine.block(model, expose = to_vb),
            init_loc_fn = init_to_value(values=init_dic)))
              
    losses = fit_with_SVI(model, guide, [x, y, n], iterations = 300, loss_tol=0.)

    # need to get per junction likelihoods. need to monte carlo here now. 
    guide_tr = pyro.poutine.trace(guide).get_trace(x, y, n)
    likelihoods = pyro.poutine.replay(model, trace = guide_tr)(x, y, n)
    medians = guide.median()
    
    return losses, guide, LeafcutterFit(
            beta = medians["beta"].clone().detach(), 
            conc = medians["conc"].clone().detach(),
            likelihoods = likelihoods,
            loss = losses[-1],
            exit_status = 0
        ) 

def svi_bbb_glm(x_null, x_full, y, n, beta_null_init, beta_full_init, conc_max = 3000., concShape=1.0001, concRate=1e-4, multiconc = True, beta_scale = np.inf, prior_prob = 0.1, eps = 1.0e-8, mc_samples = 30): 
    """
    Fits a Dirichlet binomial generalized linear model.

    Args:
        x (torch.Tensor): The input data with shape (N, P) where N is the number of observations and P is the number of predictors.
        y (torch.Tensor): The target data with shape (N, J) where J is the number of splice junctions.
        beta_init (torch.Tensor): The initial values for the coefficients, shape=(P,J)
        conc_max (float, optional): Maximum concentration parameter value. Defaults to 3000.
        concShape (float, optional): Shape parameter for concentration prior. Defaults to 1.0001.
        concRate (float, optional): Rate parameter for concentration prior. Defaults to 1e-4.
        multiconc (bool, optional): Whether to use multiple concentration parameters. Defaults to True.
        fitter (function, optional): The fitting function to use. Defaults to fit_with_lbfgs.
        guide_type (class, optional): The type of guide to use. Defaults to Simpl.
        eps (float): small pseudocount added to DM parameter for numerical stability. 

    Returns:
        LeafcutterFit: An object containing fitted model parameters, losses, and exit status.
    """
    pyro.clear_param_store()

    (N,P_null) = x_null.shape
    (N,P_full) = x_full.shape
    J = y.shape[1]
    
    model = BayesianBetaBinomialModel(P_null, P_full, J, gamma_shape = concShape, gamma_rate = concRate, multiconc = multiconc, beta_scale = beta_scale, prior_prob = prior_prob, eps = eps)
    guide = AutoGuideList(model)
    to_vb = ["conc", "beta_null", "beta_full"]
    init_dic = { # TODO: GPU 
        "beta_full": beta_full_init,
        "beta_null": beta_null_init,
        "conc": torch.full([J], 10.)
    }
    guide.add(AutoDelta(
        poutine.block(model, hide = to_vb + ['assignment']),
        init_loc_fn = init_to_value(values=init_dic)))
    guide.add(
        AutoDiagonalNormal(
            poutine.block(model, expose = to_vb),
            init_loc_fn = init_to_value(values=init_dic)))

    data = [x_null, x_full, y, n]
    losses = fit_with_SVI(
        model, guide, data, iterations = 300, loss_tol=0., 
        loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting=1))

    # need to get per junction likelihoods. need to monte carlo here now. 
    # guide_tr = pyro.poutine.trace(guide).get_trace(*data)
    # likelihoods = pyro.poutine.replay(model, trace = guide_tr)(*data)
    medians = guide.median()

    def get_sample(*nodes):
        """helper function to get discrete samples"""
        guide_trace = poutine.trace(guide).get_trace(*data)  # record the globals
        trained_model = poutine.replay(model, trace=guide_trace)  # replay the globals
        inferred_model = infer_discrete(trained_model, temperature=1, first_available_dim=-2)  # avoid conflict with data plate
        trace = poutine.trace(inferred_model).get_trace(*data)
        return([trace.nodes[node]["value"] for node in nodes] )

    samples = [get_sample('assignment') for _ in range(mc_samples)] # sample discrete RVs from approximate posterior

    stacked = [ torch.stack(g) for g in zip(*samples) ]
    kstacked = dict(zip(['assignment'], stacked))
    post_mean = {k:v.mean(0).detach().cpu().numpy() for k,v in kstacked.items()}

    return losses, guide, medians, post_mean
    
def simple_simulation(N, P, J, total_count = 100, conc = 10.):
    """
    Very simple simulation of data . 
    
    Args:
        N (int): Number of samples.
        P (int): Number of covariates.
        J (int): Number of junctions.
        total_count (int): Total count for the cluster (default: 100).
        conc (float): Concentration parameter for the Beta binomial distribution (default: 10.0).

    Returns:
        tuple: A tuple containing the following elements:
            - x (torch.Tensor): Simulated covariate data.
            - y (torch.Tensor): Simulated observed data.
            - true_beta_norm (torch.Tensor): True beta parameters after normalization.
            - g (torch.Tensor): Probabilities computed using the softmax function.

    """
    x = torch.randn((N,P))
    x[:,0] = 1. # intercept
    b = torch.randn((P,J))
    g = (x @ b).sigmoid()
    n = dist.Poisson(total_count).sample([N,J])
    bb = dist.BetaBinomial(g * conc, (1.-g)*conc, total_count = n)
    y = bb.sample()
    return(x,y,b,g)

def get_init(x, y, n, init, fitter = fit_with_lbfgs): 
    torch_types = { "device" : y.device, "dtype" : y.dtype }
    J = y.shape[1]
    if init == "brr":
        return brr_initialization(x,y,n)
    elif init == "rr": 
        return rr_initialization(x,y,n)
    elif init == "bin": # doesn't speed things up
        return fit_binomial_glm(x_null, y, n, fitter = fitter) 
    elif init == "0": 
        return torch.zeros(x.shape[1], J, **torch_types)
    else: 
        raise Exception(f"Unknown initialization strategy {init}")


def beta_binomial_anova(x_full, x_null, y, n, init = "brr", shape = None, rate = 0.5, **kwargs): 
    """
    Perform beta-binomial ANOVA analysis.

    Args:
        x_full (torch.Tensor): Full (i.e. alternative hypothesis) covariate data.
        x_null (torch.Tensor): Null (i.e. null hypothesis) covariate data.
        y (torch.Tensor): Observed junction counts
        n (torch.Tensor): ATSE total counts. 
        init (str): initialization strategy. One of "brr" (Bayesian ridge regression), "rr" (ridge regression), "bin" (binomial logistic regression) or "0" (set to 0). 
        concShape (float): Shape parameter for concentration priors (default: 1.0001).
        concRate (float): Rate parameter for concentration priors (default: 1e-4).
        multiconc (bool): Indicates whether to use separate concentration parameters for each junction (default: False).
        fitter (function): A function used to fit the model (default: fit_with_lbfgs).
        guide_type: one of BasicGuide, CleverGuide or DamCleverGuide

    Returns:
        tuple: A tuple containing the following elements:
            - loglr (float): Log-likelihood ratio statistic.
            - df (int): Degrees of freedom.
            - lrtp (float): Likelihood ratio test p-value.
            - null_fit (object): Result of fitting the null model.
            - full_fit (object): Result of fitting the full model.
            - refit_null_flag (bool): Flag indicating if the null model was refitted based on the full model.
    """
    
    (N,P_full) = x_full.shape
    (N,P_null) = x_null.shape
    J = y.shape[1]
    torch_types = { "device" : y.device, "dtype" : y.dtype }

    # fit null model
    start_time = time.time()
    beta_init_null = get_init(x_null, y, n, init)
    null_fit = fit_bb_glm(x_null, y, n, beta_init_null, **kwargs )
    elapsed_time = time.time() - start_time
    #print(f"Elapsed time: {elapsed_time} seconds, {null_fit.loss}")

    # fit full model, initialized at null model solution
    beta_init_full = torch.cat( (null_fit.beta, torch.zeros((P_full - P_null, J), **torch_types)) )
    full_fit = fit_bb_glm(x_full, y, n, beta_init_full, **kwargs )
    
    # fit full model, initialized "smartly", and check if this gives a better fit
    beta_init_full = get_init(x_full, y, n, init)
    full_fit_smart = fit_bb_glm( x_full, y, n, beta_init_full, **kwargs )
    
    smart_better = full_fit_smart.likelihoods > full_fit.likelihoods
    full_fit.beta[:,smart_better] = full_fit_smart.beta[:,smart_better]
    full_fit.likelihoods[smart_better] = full_fit_smart.likelihoods[smart_better]
    full_fit.conc[smart_better] = full_fit_smart.conc[smart_better]
    
    df = P_full - P_null
    refit_null_flag=False
    loglr = (full_fit.likelihoods - null_fit.likelihoods).detach().cpu().numpy()
    #lrtp = scipy.stats.chi2(df).sf(2.*loglr) # sf = 1-cdf
    if shape is None: 
        shape = 0.5 * df
    lrtp = scipy.stats.gamma.sf(2.*loglr, shape, scale = 1. / rate) 
    
    pretty_sig = lrtp < 0.05 # if result looks highly significant, check if we could improve fit initializing null based on full
    beta_init_null = full_fit.beta[:P_null,pretty_sig]
    refit_null = fit_bb_glm(x_null, y[:,pretty_sig], n[:,pretty_sig], beta_init_null, **kwargs )
    refit_null_flag = refit_null.likelihoods > null_fit.likelihoods[pretty_sig]
    null_fit.beta[:,pretty_sig][:,refit_null_flag] = refit_null.beta[:,refit_null_flag]
    null_fit.conc[pretty_sig][refit_null_flag] = refit_null.conc[refit_null_flag]
    
    loglr = (full_fit.likelihoods - null_fit.likelihoods).detach().cpu().numpy() 
    #lrtp = scipy.stats.chi2(df).sf(2.*loglr) # sf = 1-cdf
    lrtp = scipy.stats.gamma.sf(2.*loglr, shape, scale = 1. / rate) 
    
    return loglr, df, lrtp, null_fit, full_fit, refit_null_flag

def parametric_bootstrap(x_full, x_null, y, n, init = "brr", nboot=1, eps = 1e-6, **kwargs):
    
    J = y.shape[1]

    # fit null model
    beta_init_null = get_init(x_null, y, n, init)
    null_fit = fit_bb_glm(x_null, y, n, beta_init_null, **kwargs )

    # parametric bootstrap
    xb = x_null @ null_fit.beta
    bb = dist.BetaBinomial(
        torch.sigmoid(xb) * null_fit.conc + eps, 
        torch.sigmoid(-xb) * null_fit.conc + eps,
        total_count = n)

    boots = []
    for i in range(nboot): 
        print(i, end="")
        y_pboot = bb.sample()
        # fit on bootstrap data
        loglr_boot, _, _, _, _, _ = beta_binomial_anova(x_full, x_null, y_pboot, n, init = init, **kwargs)
        boots.append(loglr_boot)
    print("")
    
    return boots, null_fit
    #correction_factor = 2. * loglr_boot.mean() # alternatively median should be 0.454
    # correction_factor = np.median(2. * loglr_boot) / 0.454
    #m = np.mean(2. * loglr_boot)
    #v = np.var(2. * loglr_boot)
    #shape = m*m / v
    #rate = m / v

    #print(f"{shape:.3f} {rate:.3f}")
    #print(f"Parameteric bootstrap estimate of Bartlett correction {correction_factor:.3f}")

    # full fit on real data
    #return beta_binomial_anova(x_full, x_null, y, n, init = init, shape = shape, rate = rate, **kwargs)
    
    