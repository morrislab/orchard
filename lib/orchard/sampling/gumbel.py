########################################################################
# gumbel.py
#
# Contains the source code for performing the Gumbel-Max trick
# and computing 'total perturbed log-likelihoods'.
########################################################################
import numpy as np

def gumbel(generator, shape=1):
    """Samples a Gumbel random variable
    
    Parameters
    ----------
    generator : object
        a numpy default_rng object used for reproducible random sampling
    shape : int, optional
        the number of Gumbel random variables to sample (shape>1 returns an array)

    Returns
    -------
    float or ndarray
        depending on the shape parameter, the function will either return a single float or an array of floats
    """
    return _gumbel(generator.random(shape))

def _gumbel(u):
    """Computes a Gumbel random variable
    
    Parameters
    ----------
    u : float or ndarray
        the value(s) to transform into a Gumbel random variable

    Return
    ------
    float or ndarray
        depending on the size of the input u, the function will either return a single float or an array of floats
    """
    return -np.log(-np.log(u))

def gumbel_with_maximum(phi, T, generator):
    """Shifts Gumbel ranodom variable(s) by the maximum
    
    Parameters
    ----------
    phi : ndarray
        an array of log-likelihoods
    T : float
        previous total perturbed log-likelihood
    generator : object
        a numpy default_rng object used for reproducible random sampling

    Returns
    -------
    ndarray
        an array of total perturbed log-likelihoods for each original value of phi (log-likelihoods)        
    """
    g_phi = phi + gumbel(generator, phi.shape)
    Z = np.max(g_phi)
    return _shift_by_maximum(g_phi, T, Z)

def _shift_by_maximum(g_phi, T, Z):
    """Shifts the Gumbel random variables conditionally on their maximums, and computes the 
    total perturbed log-likelihood for each input in g_phi

    Parameters
    ----------
    g_phi : ndarray
        an array of Gumbel random variables (i.e., perturbed log-likelihoods)
    T : float
        previous total perturbed log-likelihood
    Z : float
        the maximum on which all g_phi are conditional on

    Returns
    -------
    ndarray
        an array of total perturbed log-likelihoods for each input in g_phi
    """
    g_hat = 1-np.exp(g_phi - Z)

    # Watch out here too! np.log(g_hat) may return -inf
    # and if you instead do something like np.log(np.maximum(1e-300, g_hat))
    # you might end up discarding the solution with the best log probability!
    u = T - np.add(g_phi, np.log(g_hat, where=g_hat>0))

    # BEWARE, we can get some numerical stability issues
    # For example if we compute u as follows:
    #   u = T - g_phi + np.log(g_hat, where=g_hat!=0))
    # 
    # Python addition will take place between g_phi and the 
    # log of g_hat, and this can lead to nan values
    # 
    # By switching to use np.add, we ensure the types of the
    # two arrays are the same and avoid any numerical instability!
    if any(np.isnan(u)):
        raise Exception("Numerical stability issue found when computing gumbel variable")
        
    return T - np.maximum(0, u) - np.log1p(np.exp(-np.abs(u)))
