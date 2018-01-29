__author__ = 'Haohan Wang'

import scipy.linalg as linalg
import scipy
import numpy as np
from scipy import stats

def matrixMult(A, B):
    try:
        linalg.blas
    except AttributeError:
        return np.dot(A, B)

    if not A.flags['F_CONTIGUOUS']:
        AA = A.T
        transA = True
    else:
        AA = A
        transA = False

    if not B.flags['F_CONTIGUOUS']:
        BB = B.T
        transB = True
    else:
        BB = B
        transB = False

    return linalg.blas.dgemm(alpha=1., a=AA, b=BB, trans_a=transA, trans_b=transB)

def factor(X, rho):
    """
    computes cholesky factorization of the kernel K = 1/rho*XX^T + I

    Input:
    X design matrix: n_s x n_f (we assume n_s << n_f)
    rho: regularizaer

    Output:
    L  lower triangular matrix
    U  upper triangular matrix
    """
    n_s, n_f = X.shape
    K = 1 / rho * scipy.dot(X, X.T) + scipy.eye(n_s)
    U = linalg.cholesky(K)
    return U

def tstat(beta, var, sigma, q, N, log=False):

    """
       Calculates a t-statistic and associated p-value given the estimate of beta and its standard error.
       This is actually an F-test, but when only one hypothesis is being performed, it reduces to a t-test.
    """
    ts = beta / np.sqrt(var * sigma)
    # ts = beta / np.sqrt(sigma)
    # ps = 2.0*(1.0 - stats.t.cdf(np.abs(ts), self.N-q))
    # sf == survival function - this is more accurate -- could also use logsf if the precision is not good enough
    if log:
        ps = 2.0 + (stats.t.logsf(np.abs(ts), N - q))
    else:
        ps = 2.0 * (stats.t.sf(np.abs(ts), N - q))
    if not len(ts) == 1 or not len(ps) == 1:
        raise Exception("Something bad happened :(")
        # return ts, ps
    return ts.sum(), ps.sum()

def nLLeval(ldelta, Uy, S, REML=True):
    """
    evaluate the negative log likelihood of a random effects model:
    nLL = 1/2(n_s*log(2pi) + logdet(K) + 1/ss * y^T(K + deltaI)^{-1}y,
    where K = USU^T.

    Uy: transformed outcome: n_s x 1
    S:  eigenvectors of K: n_s
    ldelta: log-transformed ratio sigma_gg/sigma_ee
    """
    n_s = Uy.shape[0]
    delta = scipy.exp(ldelta)

    # evaluate log determinant
    Sd = S + delta
    ldet = scipy.sum(scipy.log(Sd))

    # evaluate the variance
    Sdi = 1.0 / Sd
    Uy = Uy.flatten()
    ss = 1. / n_s * (Uy * Uy * Sdi).sum()

    # evalue the negative log likelihood
    nLL = 0.5 * (n_s * scipy.log(2.0 * scipy.pi) + ldet + n_s + n_s * scipy.log(ss))

    if REML:
        pass

    return nLL

def testRepresentability(data):
    l = []
    [m, n] = data.shape
    for i in range(n):
        X1 = data[:, i]
        X2 = np.delete(data, i, 1)
        C11 = np.dot(X1.T, X1) * 1.0 / n
        C21 = np.dot(X2.T, X1) * 1.0 / n
        if C11 == 0:
            c = 0
        else:
            ii = 1.0 / C11
            r = np.abs(np.dot(C21, ii))
            c = len(np.where(r >= 1+1e-5)[0])
        l.append(c)
    return np.array(l)

def solve_cg(A, b, x, args, tol=10**(-8), k_max=None):
    if k_max == None:
        k_max = x.shape[0]
    k = 0
    r = b - A(x, **args)
    rho_0 = np.dot(r, r)
    rho_1 = rho_0
    while (rho_1 > tol) and (k < k_max):
        k += 1
        if k == 1:
            p = r
        else:
            beta = rho_1 / rho_0
            p = r + beta * p
        w = A(p, **args)
        alpha = rho_1 / np.dot(p, w)
        x = x + alpha * p
        r = r - alpha * w
        rho_0 = rho_1
        rho_1 = np.dot(r, r)
    return x, k

def A_trace(x, M, D):
    return np.dot(M.T, np.dot(M, x)) + D * x

def norm(w):
    return np.sqrt(np.sum(w**2))

# def logisticRegressionGradientSolver(w, X, xx, xy, D, lr, tol, maxIter):
#     resi_prev = np.inf
#     xxd = xx + D.T
#     tmp = 1/(1+np.exp(-np.dot(X, w)))
#     diff = xy-np.dot(X.T, tmp)-np.dot(D, w)
#     resi = logisticRegressionSolverCost(diff)
#     step = 0
#     while np.abs(resi_prev - resi) > tol and step < maxIter:
#         keepRunning = True
#         resi_prev = resi
#         runningStep = 0
#         while keepRunning and runningStep < 10:
#             runningStep += 1
#             prev_w = w
#             grad = -np.dot(xxd, diff)
#             w = w - grad * lr
#             tmp = 1/(1+np.exp(-np.dot(X, w)))
#             diff = xy-np.dot(X.T, tmp)-np.dot(D, w)
#             keepRunning = stopCheck(prev_w, w, grad, xxd, xy)
#             # keepRunning = False
#             if keepRunning:
#                 lr = 0.5 * lr
#
#         step += 1
#         resi = logisticRegressionSolverCost(diff)
#         if resi > resi_prev:
#             return np.zeros_like(w)
#     return w
#
def logisticRegressionSolverCost(t):
    return np.linalg.norm(t)

def logisticGradient(X, D, Xi, diff):
    return np.dot(X.transpose(), diff) - np.multiply(D, np.dot(diff, Xi.T))

# def stopCheck(prev, new, pg, X, y):
#     if np.linalg.norm(y - np.dot(X, new)) <= \
#                             np.linalg.norm(y - np.dot(X, new)) + np.sum(np.dot(pg.transpose(), (new - prev))):
#         return False
#     else:
#         return True

def logisticRegressionGradientSolver(w, X, y, D, lr, tol, maxIter, quiet=True):
    resi_prev = np.inf
    xw = np.dot(X, w).T
    xw = np.minimum(xw, 500)
    tmp = 1/(1+np.exp(xw))
    xi = np.dot(X.T, np.linalg.pinv(np.dot(X, X.T)))

    Dxw = np.dot(np.multiply(D, w), xi)
    diff = y*xw - tmp - Dxw
    resi = logisticRegressionSolverCost(diff)

    step = 0
    while resi_prev > resi and resi_prev - resi > tol and step < maxIter:
        resi_prev = resi
        if not quiet:
            print resi_prev
        # keepRunning = True
        # runningStep = 0
        # while keepRunning and runningStep < 10:
        #     runningStep += 1
        prev_w = w
        grad = logisticGradient(X, D, xi, diff)
        w = w - grad * lr
        xw = np.dot(X, w).T
        xw = np.minimum(xw, 500)
        tmp = 1/(1+np.exp(xw))
        Dxw = np.dot(np.multiply(D, w), xi)
        diff = y - tmp - Dxw
            # keepRunning = stopCheck(prev_w, w, grad, X, y)
            # keepRunning = False
            # if keepRunning:
            #     lr = 0.5 * lr

        step += 1
        resi = logisticRegressionSolverCost(diff)
        if resi > resi_prev:
            return w
    return w