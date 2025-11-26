import torch
import pyro
import pyro.distributions as dist
from pyro.infer.autoguide import AutoDiagonalNormal, AutoGuideList, AutoDelta, init_to_value, AutoDiscreteParallel
from pyro.infer import SVI, Trace_ELBO, config_enumerate, infer_discrete
from pyro.ops.indexing import Vindex
from pyro import poutine

import torch
import pyro
import pyro.distributions as dist
from pyro.infer import config_enumerate
from tqdm import tqdm
from sklearn import linear_model
import numpy as np

from leafcutter.differential_splicing.bb_glm import brr_initialization

pyro.enable_validation(False) # for binomial mixture https://github.com/pyro-ppl/pyro/issues/3419

def convertr(hyperparam, name, device = "cpu"): 
    #return pyro.sample(name, hyperparam) if (type(hyperparam) != float) else torch.tensor(hyperparam, device = device)
    is_dist = isinstance(hyperparam, pyro.distributions.Distribution)
    return pyro.sample(name, hyperparam) if (is_dist) else (
        torch.tensor(hyperparam, device = device) if (type(hyperparam) == float) else hyperparam
    )

def simulate_data(
    N = 6, # samples
    P = 3, # covariates
    J = 200, # junctions
    conc = 10,
    prop_non_full = 0.5
): 
    x_null = torch.randn([N,P-1])
    half_N = N // 2
    column = torch.cat([torch.zeros(half_N), torch.ones(N - half_N)]).unsqueeze(1)
    x_full = torch.cat([x_null, column], dim=1)
    
    b = torch.rand([P,J]).sign() * (torch.randn((P,J)) + 1.)
    g_full = (x_full @ b).sigmoid()
    g_null = (x_null @ b[:-1,:]).sigmoid()
    is_full = torch.rand(J) < prop_non_full
    b[-1,~is_full] = 0. 
    #g = (x_full @ b).sigmoid()
    g = torch.where(is_full, g_full, g_null)
    n = dist.Poisson(100).sample([N,J])
    y = dist.BetaBinomial(g * conc, (1.-g)*conc, total_count = n).sample()

    return x_null, x_full, y, n, is_full, b

class BayesianBetaBinomialModel(pyro.nn.PyroModule):
    
    def __init__(self, eps = 1e-8, gamma_shape = 2., gamma_rate = 0.2, beta_scale = 1., multiconc = True): 
        """
        Initialize the BayesianBetaBinomialModel.

        For gamma_shape, gamma_rate, beta_scale, and prior_prob if the parameter is a distribution object, then it will be 
        learned under that prior. 

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
        assert(multiconc)
        self.eps = eps
        self.binomial = gamma_shape is None
        self.gamma_shape = gamma_shape
        self.gamma_rate = gamma_rate
        self.multiconc = multiconc
        self.beta_scale = beta_scale # maybe this should be learned per covariate? 

    def forward(self, x, y, n):

        P = x.shape[1]
        N, J = y.shape
        device = x.device

        beta_scale = convertr(self.beta_scale, "beta_scale", device = device)

        if not self.binomial:
            gamma_shape = convertr(self.gamma_shape, "gamma_shape", device = device)
            gamma_rate = convertr(self.gamma_rate, "gamma_rate", device = device)

        with pyro.plate("P", P):
            beta = pyro.sample("beta", dist.Laplace(0., beta_scale).expand([P, J]).to_event(1))

        logits = x @ beta # N x J

        with pyro.plate("J", J):

            if self.binomial:
                bb = dist.Binomial(logits = logits.T, total_count = n.T).to_event(1)
            else:
                conc = pyro.sample("conc", dist.Gamma(gamma_shape, gamma_rate)) # J
                g = logits.sigmoid()
                bb = dist.BetaBinomial(
                    concentration1=(g * conc + self.eps).T,
                    concentration0=((1 - g) * conc + self.eps).T,
                    total_count=n.T
                ).to_event(1)

            with pyro.plate("N", N):
                pyro.sample("obs", bb, obs=y.T)
            
    def fit(self, x, y, n, beta_init = None, lr=0.01, iterations = 500):

        N, J = y.shape
        pyro.clear_param_store()

        if beta_init is None:
            beta_init = brr_initialization(x, y, n)

        init_dic = { # TODO: GPU 
            "beta": beta_init,
            "conc": torch.full([J], 10.), 
            "gamma_shape" : 2., 
            "gamma_rate" : 0.2, 
            "beta_scale" : torch.std(beta_init)
        }

        # attempt to calculate ELBO
        guide = AutoDiagonalNormal(self, init_loc_fn = init_to_value(values=init_dic))
        elbo_func = pyro.infer.Trace_ELBO()
        optim = pyro.optim.Adam({"lr": lr})
        svi = SVI(self, guide, optim, loss=elbo_func)

        # Optimization
        losses = []
        for step in tqdm(range(iterations)):
            loss = svi.step(x, y, n)
            losses.append(loss)

        return losses, guide

def bin_then_bb_glm(x, y, n, lr = 0.01, iterations = 500, gamma_shape = 2., num_particles = 30, **kwargs):

    print("Initial binomial GLM fitting.")
    binomial_glm = BayesianBetaBinomialModel(gamma_shape = None, **kwargs)
    losses, guide = binomial_glm.fit(x, y, n, lr=lr, iterations=iterations) # fit binomial glm first, seems worth it here
    beta_init = guide.median()['beta']

    print("Beta-binomial GLM fitting.")
    bb_glm = BayesianBetaBinomialModel(gamma_shape = gamma_shape, **kwargs)
    bb_losses, guide = bb_glm.fit(x, y, n, beta_init = beta_init, iterations=iterations, **kwargs) # initialize beta binomial glm from there
    final_elbo = pyro.infer.Trace_ELBO(num_particles = num_particles)(bb_glm, guide)(x, y, n).item()
    return losses+bb_losses, final_elbo, guide.median()["beta"], guide.median()["conc"]

def estimate_marginal_posterior(logw, alpha, pi = None): 

    num_samples = logw.shape[0]
    if alpha == 1.: # variational inference
        log_marg = logw.mean(0)
    elif alpha == 0.: # importance sampling
        log_marg = logw.logsumexp(0) - torch.log(torch.tensor(num_samples))
    else: # renyi
        one_minus_alpha = 1. - alpha
        log_marg = (one_minus_alpha * logw).logsumexp(0) - torch.log(torch.tensor(num_samples))
        log_marg /= one_minus_alpha
    log_bayes_factor = log_marg[1] - log_marg[0]
    if pi is not None:
        log_marg += pi
    log_marg -= log_marg.logsumexp(0, keepdim = True)
    prob = log_marg.exp()
    return prob[1], log_bayes_factor

def effective_sample_size(logw):
    # \text{ESS} = \frac{(\sum_{i=1}^N w_i)^2}{\sum_{i=1}^N w_i^2}
    log_numerator = 2. * logw.logsumexp(0)
    log_denominator = (2.*logw).logsumexp(0)
    return torch.exp(log_numerator - log_denominator)
        
class SpikeAndSlabModel(pyro.nn.PyroModule):
    
    def __init__(
            self, 
            eps = 1e-8, 
            gamma_shape = 2., 
            gamma_rate = .2, 
            beta_scale = 1., 
            prior_prob = torch.tensor([0.9,0.1]), 
            multiconc = True, 
            per_hyp_conc = True): 
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
            prior_prob (float): Probability vector P(null), P(full) or K=2 Dirichlet prior
            per_hyp_conc (bool): Indicates whether to use a separate concentration parameter for each hypothesis.
        """
        super().__init__()
        assert(multiconc) # False not implemented
        self.eps = eps
        self.binomial = gamma_shape is None
        self.gamma_shape = gamma_shape
        self.gamma_rate = gamma_rate
        self.multiconc = multiconc
        self.beta_full_scale = beta_scale
        self.beta_null_scale = beta_scale
        self.prior_prob = prior_prob
        self.per_hyp_conc = per_hyp_conc
        self.guide = None

    def forward(self, x_null, x_full, y, n):

        P_full = x_full.shape[1]
        P_null = x_null.shape[1]
        N, J = y.shape
        device = x_full.device

        beta_full_scale = convertr(self.beta_full_scale, "beta_full_scale", device = device)
        beta_null_scale = convertr(self.beta_null_scale, "beta_null_scale", device = device)

        with pyro.plate("P_null", P_null):
            beta_null = pyro.sample("beta_null", dist.Laplace(0., beta_null_scale).expand([P_null, J]).to_event(1))

        with pyro.plate("P_full", P_full):
            beta_full = pyro.sample("beta_full", dist.Laplace(0., beta_full_scale).expand([P_full, J]).to_event(1))

        # should these be different for null vs full? maybe not since should capture technical only
        gamma_shape = convertr(self.gamma_shape, "gamma_shape", device = device)
        gamma_rate = convertr(self.gamma_rate, "gamma_rate", device = device)

        prior_prob = convertr(self.prior_prob, "prior_prob")
        mix = dist.Categorical(prior_prob).expand([J])
        #print("mix shapes", mix.batch_shape, mix.event_shape) # [J],[]

        logits_null = x_null @ beta_null
        logits_full = x_full @ beta_full

        logits_combined = torch.stack((logits_null, logits_full), dim=0)  # 2 x N x J
        #print("logits_combined.shape", logits_combined.shape)

        g = logits_combined.sigmoid()

        with pyro.plate("J", J):
            conc_dist = dist.Gamma(gamma_shape, gamma_rate)
            if self.per_hyp_conc:
                conc_param = pyro.sample("conc", conc_dist.expand([2,J]))[:,None,:] # 2 x J
            else: 
                conc_param = pyro.sample("conc", conc_dist.expand([J])) # J
            #print(conc_param.shape)

            bb = dist.BetaBinomial(
                concentration1=(g * conc_param + self.eps).permute(2, 0, 1), # 2 x N x J -> J x 2 x N
                concentration0=((1 - g) * conc_param + self.eps).permute(2, 0, 1),
                total_count=n.T[:,None,:] # N x J -> J x 1 x N so broadcasts over components
            )
            #print("bb shapes:",bb.batch_shape, bb.event_shape) # [J x K=2 x N], []
            comp = dist.Independent(bb, reinterpreted_batch_ndims = 1)
            #print("comp shapes:", comp.batch_shape, comp.event_shape) # [J x K=2], [N]

            mixture = dist.MixtureSameFamily(
                mix,
                comp
            ) # .to_event(1) isn't needed because of the J plate
            #print("mixture shapes:", mixture.batch_shape, mixture.event_shape) # [], [J x N]
            with pyro.plate("N", N):
                pyro.sample("obs", mixture, obs=y.T)

    def fit(
        self, 
        x_null,
        x_full, 
        y, 
        n, 
        beta_null_init = None, 
        beta_full_init = None, 
        conc_null_init = None, 
        conc_full_init = None,
        alpha = 1., 
        num_particles = 1, 
        lr=0.01, 
        iterations = 1000
        ):

        pyro.clear_param_store()

        losses_null = None
        if beta_null_init is None:
            print("Fitting null model")
            losses_null, final_elbo, beta_null_init, conc_null_init = bin_then_bb_glm(x_null, y, n, lr = lr, iterations = iterations // 2) 

        losses_full = None
        if beta_full_init is None:
            print("Fitting full model")
            losses_full, final_elbo, beta_full_init, conc_full_init = bin_then_bb_glm(x_full, y, n, lr = lr, iterations = iterations // 2)

        if conc_null_init is None:
            conc_null_init = torch.full([J], 10.)

        if conc_full_init is None:
            conc_full_init = torch.full([J], 10.)

        if self.per_hyp_conc:
            conc_init = torch.stack([conc_null_init, conc_full_init])
        else: 
            conc_init = conc_null_init

        init_dic = { # TODO: GPU 
            "beta_full": beta_full_init,
            "beta_null": beta_null_init,
            "conc": conc_init, 
            "mixing_probs" : torch.tensor([0.9,0.1])
        }

        guide = AutoGuideList(self)
        guide.add(AutoDiagonalNormal(
            poutine.block(self, expose = ['beta_null']),
            init_loc_fn = init_to_value(values=init_dic)))
        guide.add(AutoDiagonalNormal(
            poutine.block(self, expose = ['beta_full']),
            init_loc_fn = init_to_value(values=init_dic)))
        guide.add(AutoDiagonalNormal(
            poutine.block(self, expose = ['conc']),
            init_loc_fn = init_to_value(values=init_dic)))
        if isinstance(self.prior_prob, torch.distributions.Distribution):
            guide.add(AutoDiagonalNormal(
                poutine.block(self, expose = ['mixing_probs']),
                init_loc_fn = init_to_value(values=init_dic)))
        
        if isinstance(self.beta_full_scale, torch.distributions.Distribution):
            guide.add(AutoDiagonalNormal( 
                poutine.block(self, expose = ['beta_full_scale']),
                init_loc_fn = init_to_value(values=init_dic)))
        if isinstance(self.beta_null_scale, torch.distributions.Distribution):
            guide.add(AutoDiagonalNormal( 
                poutine.block(self, expose = ['beta_null_scale']),
                init_loc_fn = init_to_value(values=init_dic)))
        if isinstance(self.gamma_shape, torch.distributions.Distribution):
            guide.add(AutoDiagonalNormal( 
                poutine.block(self, expose = ['gamma_shape']),
                init_loc_fn = init_to_value(values=init_dic)))
        if isinstance(self.gamma_rate, torch.distributions.Distribution):
            guide.add(AutoDiagonalNormal( 
                poutine.block(self, expose = ['gamma_rate']),
                init_loc_fn = init_to_value(values=init_dic)))

        self.guide = guide # AutoDiagonalNormal(self, init_loc_fn = init_to_value(values=init_dic))

        loss_func = pyro.infer.Trace_ELBO(num_particles = num_particles) if (
            alpha == 1.) else pyro.infer.RenyiELBO(alpha=alpha, num_particles=num_particles) 

        optim = pyro.optim.Adam({"lr": lr})
        svi = SVI(self, self.guide, optim, loss=loss_func)

        print("Fitting spike and slab joint model")
        # Optimization
        losses = []
        for step in tqdm(range(iterations)):
            loss = svi.step(x_null, x_full, y, n)
            losses.append(loss)

        return losses_null, losses_full, losses

    def get_posterior_map(self, x_null, x_full, y, n): 

        N, J = y.shape

        guide_median = self.guide.median()
        beta_full = guide_median["beta_full"]
        beta_null = guide_median["beta_null"]
        conc_param = guide_median["conc"]
        
        mixing_probs = torch.tensor([0.9,0.1]) # TODO: learned
        mix = dist.Categorical(mixing_probs).expand([J])
        
        logits_full = x_full @ beta_full
        logits_null = x_null @ beta_null
        logits_combined = torch.stack((logits_null, logits_full), dim=0)  # 2 x N x J
        
        g = logits_combined.sigmoid()

        if self.per_hyp_conc:
            conc_param = conc_param[:,None,:]
        bb = dist.BetaBinomial(
            concentration1=(g * conc_param + 1e-8).permute(2, 0, 1), # 2 x N x J
            concentration0=((1 - g) * conc_param + 1e-8).permute(2, 0, 1),
            total_count=n.T[:,None,:] # N x J -> J x 1 x N so broadcasts over components
        )
        
        comp = dist.Independent(bb, reinterpreted_batch_ndims = 1)

        mixture = dist.MixtureSameFamily(mix, comp) # .to_event(1) this would get us just one log_prob
        
        log_prob_x = mixture.component_distribution.log_prob(mixture._pad(y.T))  # [S, B, k]
        log_mix_prob = torch.log_softmax(mixture.mixture_distribution.logits, dim=-1)  # [B, k]
        log_prob = log_prob_x + log_mix_prob
        log_prob -= log_prob.logsumexp(-1, keepdim=True)
        probs = log_prob.exp()

        return probs[:,1].detach()

    def get_importance_sampling_weights(
        self,
        x_null,
        x_full, 
        y, 
        n,
        num_samples = 100,
        sample_conc = True
    ):
        P_full = x_full.shape[1]
        P_null = x_null.shape[1]
        N, J = y.shape
        
        conc = self.guide.median()["conc"]
        
        logw = []
        with torch.no_grad(): 
            for _ in range(num_samples):
                guide_tr = poutine.trace(self.guide).get_trace()
                conditioned_model = poutine.replay(lambda: self(x_null, x_full, y, n), guide_tr)
                model_tr = poutine.trace(conditioned_model).get_trace()
            
                model_tr.nodes["obs"]["fn"].log_prob(model_tr.nodes["obs"]["value"]) # N x J
            
                beta_null = guide_tr.nodes["beta_null"]["value"]
                beta_full = guide_tr.nodes["beta_full"]["value"]
                
                mixture = model_tr.nodes["obs"]["fn"]
                logp_obs = mixture.component_distribution.log_prob(mixture._pad(y.T))  # [S, B, k]
                #logp_mix = torch.log_softmax(mixture.mixture_distribution.logits, dim=-1)  # do we want to include this? 
                #log_prob = logp_obs + logp_mix # N x J x 2
            
                logp_beta_full = model_tr.nodes["beta_full"]["fn"].base_dist.log_prob(beta_full).sum(0)
                logp_beta_null = model_tr.nodes["beta_null"]["fn"].base_dist.log_prob(beta_null).sum(0) # 2 x J
                logp_beta = torch.stack([logp_beta_null,logp_beta_full]) # 2 x J
            
                logp = logp_obs.sum(0).T + logp_beta 
            
                q_beta_null = dist.Normal( 
                    guide_tr.nodes['guide.0.loc']["value"].reshape([P_null, J]),
                    guide_tr.nodes['guide.0.scale']["value"].reshape([P_null, J])
                )
                logq_beta_null = q_beta_null.log_prob(beta_null).sum(0) # P_null x J 
            
                assert(torch.allclose(guide_tr.nodes['_guide.1_latent']["value"].reshape([P_full, J]), beta_full))
                q_beta_full = dist.Normal( 
                    guide_tr.nodes['guide.1.loc']["value"].reshape([P_full, J]),
                    guide_tr.nodes['guide.1.scale']["value"].reshape([P_full, J])
                )
                logq_beta_full = q_beta_full.log_prob(beta_full).sum(0) # P_full x J 
                logq_beta = torch.stack([logq_beta_null,logq_beta_full]) # 2 x J

                conc_shape = [2,J] if self.per_hyp_conc else [J]
                # this is log-Normal so a bit more complex H[y] = H[x] + E[log|g'(x)'] where y=g(x)). Not handled properly? 
                q_conc = dist.Normal(
                    guide_tr.nodes['guide.2.loc']["value"].reshape(conc_shape),
                    guide_tr.nodes['guide.2.scale']["value"].reshape(conc_shape)
                )
                
                logq = logq_beta
            
                if sample_conc: 
                    log_conc = guide_tr.nodes['_guide.2_latent']["value"].reshape(conc_shape)
                    conc = guide_tr.nodes["conc"]["value"]
                    logq += q_conc.log_prob(log_conc) - log_conc # explicit Jacobian!?
                    logp += model_tr.nodes["conc"]["fn"].log_prob(conc)
            
                logw.append(logp - logq)
            logw = torch.stack(logw)
            pi = torch.log_softmax(mixture.mixture_distribution.logits, dim=-1)[0].T
        return logw, pi

    def estimate_marginal_posterior(
        self, 
        x_null,
        x_full, 
        y, 
        n,
        alpha = 0., # weirdly alpha=1. seems to be the best calibrated
        num_samples = 100
    ): 
        logw, pi = self.get_importance_sampling_weights(x_null, x_full, y, n, num_samples = num_samples)
        return estimate_marginal_posterior(logw, alpha, pi = pi)
        
