import scipy.optimize as opt
from sklearn.preprocessing import normalize

from helpingMethods import *

class CMM:
    def __init__(self, lam=1.0, rho=1.0, lr=1.0, maxIter=100, tol=1e-3, maxADMMIter=100, maxPGDIter=100, logistic=False,
                 quiet=True):
        self.lam1 = lam
        self.lam2 = lam
        self.rho = rho
        self.lr = lr
        self.maxIter = maxIter
        self.tol = tol
        self.maxADMMIter = maxADMMIter
        self.maxPGDIter = maxPGDIter
        self.logistic = logistic
        self.decay = 0.5
        self.convergenceCheckerStart = 2
        self.quiet = quiet

    def setLambda(self, lam):
        self.lam1 = lam
        self.lam2 = lam

    def setLambda1(self, lam):
        self.lam1 = lam

    def setLambda2(self, lam):
        self.lam2 = lam

    def setRho(self, rho):
        self.rho = rho

    def setLearningRate(self, lr):
        self.lr = lr

    def checkConvergence(self):
        '''
        :return: True for keep running, false for stop
        '''
        if self.convergenceCheckerCount < self.convergenceCheckerStart:
            self.b1Prev = self.b1
            self.b2Prev = self.b2
            self.convergenceCheckerCount += 1
            return True
        if np.linalg.norm(self.b1 - self.b1Prev) < self.tol and np.linalg.norm(self.b2 - self.b2Prev) < self.tol:
            return False
        self.b1Prev = self.b1
        self.b2Prev = self.b2
        return True

    def rescale(self, a):
        return a / np.max(np.abs(a))

    def selectValues(self, Kva):
        r = np.zeros_like(Kva)
        n = r.shape[0]
        tmp = self.rescale(Kva[:-1])
        ind = 0
        for i in range(n / 2, n - 2):
            if tmp[i + 1] - tmp[i] > 1.0 / n:
                ind = i + 1
                break
        r[ind:] = Kva[ind:]
        r[n - 1] = Kva[n - 1]
        return r

    def estimating_variance_helper(self, y, X, S=None, U=None, numintervals=500, ldeltamin=-5, ldeltamax=5, scale=0):
        ldeltamin += scale
        ldeltamax += scale

        y = y - np.mean(y)

        # if S is None or U is None:
        S = None
        U = None
        K = matrixMult(X, X.T)
        S, U = linalg.eigh(K)
        Uy = scipy.dot(U.T, y)

        ns = X.shape[0]

        S = self.selectValues(S)
        S = normalize(S.reshape([1, ns])).reshape([ns])

        nllgrid = scipy.ones(numintervals + 1) * scipy.inf
        ldeltagrid = scipy.arange(numintervals + 1) / (numintervals * 1.0) * (ldeltamax - ldeltamin) + ldeltamin
        for i in scipy.arange(numintervals + 1):
            nllgrid[i] = nLLeval(ldeltagrid[i], Uy, S)  # the method is in helpingMethods

        nllmin = nllgrid.min()
        ldeltaopt_glob = ldeltagrid[nllgrid.argmin()]

        for i in scipy.arange(numintervals - 1) + 1:
            if (nllgrid[i] < nllgrid[i - 1] and nllgrid[i] < nllgrid[i + 1]):
                ldeltaopt, nllopt, iter, funcalls = opt.brent(nLLeval, (Uy, S),
                                                              (ldeltagrid[i - 1], ldeltagrid[i], ldeltagrid[i + 1]),
                                                              full_output=True)
                if nllopt < nllmin:
                    nllmin = nllopt
                    ldeltaopt_glob = ldeltaopt

        delta0 = scipy.exp(ldeltaopt_glob)
        Sdi = 1. / (S + delta0)
        Sdi = normalize(Sdi.reshape([1, ns])).reshape(ns)
        SUy = np.square(scipy.dot(U.T, y))
        SUy = normalize(SUy.reshape(1, ns)).reshape(ns)
        ratio = SUy * Sdi
        sigmaU = np.sum(ratio)
        sigmaE = sigmaU * delta0

        Sdi_sqrt = scipy.sqrt(Sdi)
        SUX = scipy.dot(U.T, X)
        X = SUX * scipy.tile(Sdi_sqrt, (X.shape[1], 1)).T
        X = normalize(X, axis=0)
        SUy = scipy.dot(U.T, y)
        y = SUy * scipy.reshape(Sdi_sqrt, (ns))
        y = normalize(y.reshape(1, ns)).reshape(ns)

        return np.sum(S) * sigmaU, sigmaE, X, y

    def estimatingVariance(self, S1, U1, S2, U2):
        '''
        :return: tr(K_1*sigma_u1^2), tr(K_2*sigma_u2^2), sigma_e1^2, sigma_e2^2
        '''
        tu1, se1, self.X1, self.y1 = self.estimating_variance_helper(self.y1, self.X1, S=S1, U=U1)
        tu2, se2, self.X2, self.y2 = self.estimating_variance_helper(self.y2, self.X2, S=S2, U=U2)
        return tu1, tu2, se1, se2

    def calculateSS(self, X, y, b, tu, se):
        tmp = np.dot(X, b)
        return np.dot(y.T, y) * (1.0 / y.shape[0]) + np.dot(tmp.T, tmp) * (1.0 / X.shape[0]) + 2 * (tu + se) * (
        1.0 / (y.shape[0] + X.shape[0]))

    def calculatingSigmaT(self, tu1, tu2, se1, se2):
        '''
        :return: sigma11, sigma22, t
        '''
        s11 = self.calculateSS(self.X2, self.y1, self.b1, tu1, se1)
        s22 = self.calculateSS(self.X1, self.y2, self.b2, tu2, se2)
        # t = s22 * np.linalg.norm(self.y1 - np.dot(self.X1, self.b1), ord=2)*(1.0/self.X1.shape[0]) + s11 * np.linalg.norm(
        #     self.y2 - np.dot(self.X2, self.b2), ord=2)*(1.0/self.X2.shape[0])
        s12 = np.dot(self.y1.T, np.dot(self.X1, self.b2)) * (1.0 / self.y1.shape[0]) + \
              np.dot(self.y2.T, np.dot(self.X2, self.b1)) * (1.0 / self.y2.shape[0]) + \
              (tu1 + se1 + tu2 + se2) * (1.0 / (self.y1.shape[0] + self.y2.shape[0]))
        s21 = s12
        t = s11 * s22 - s12 * s21
        t = max(t, 1e-5)
        # print t
        return s11, s22, t

    def solveBeta(self, X, y, b, b_, L, c, sign, lam):
        '''
        :return: updated beta
        '''
        self.bpg = b
        self.bpg2 = b_
        lr = self.lr
        resi_prev = np.inf
        resi = self.cost(X, y, c, L, lam)
        step = 0
        while resi_prev - resi > self.tol and step < self.maxPGDIter:
            if not self.quiet:
                print '\t\t\t\tPGD Iteration', step, resi

            resi_prev = resi
            pg = self.proximal_gradient(X, y, c, L, sign)
            self.bpg = self.proximal_proj(self.bpg - pg * lr, lr, lam)
            step += 1
            resi = self.cost(X, y, c, L, lam)
        return self.bpg

    def cost(self, X, y, c, L, lam):
        if self.logistic:
            v = (np.dot(X, self.bpg)).T
            tmp = - c * np.sum(y * v - np.log(1 + np.exp(v)))
        else:
            tmp = c * np.sum(np.square(y - np.dot(X, self.bpg)).transpose())
        tmp = tmp * (1.0 / X.shape[0]) + lam * linalg.norm(self.bpg, ord=1)
        tmp += self.rho * linalg.norm(self.bpg - self.bpg2, ord=2) + np.dot(L.T, self.bpg - self.bpg2)
        return tmp

    def proximal_gradient(self, X, y, c, L, sign):
        if self.logistic:  # this is the correct derivation of log likelihood https://beckernick.github.io/logistic-regression-from-scratch/
            tmp = - c * np.dot(X.transpose(), (y.reshape((y.shape[0], 1)) - 1. / (1 + np.exp(-np.dot(X, self.bpg)))))
        else:
            tmp = -c * np.dot(X.transpose(), (y.reshape((y.shape[0], 1)) - (np.dot(X, self.bpg))))
        tmp = tmp * (1.0 / X.shape[0])
        tmp += 2 * self.rho * (self.bpg - self.bpg2) * sign + L
        return tmp

    def proximal_proj(self, B, lr, lam):
        t = lam * lr
        zer = np.zeros_like(B)
        result = np.maximum(zer, B - t) - np.maximum(zer, -B - t)
        return result

    def stopCheck(self, prev, new, pg, X, y, lam):
        if np.square(linalg.norm((y - (np.dot(X, new))))) <= \
                                np.square(linalg.norm((y - (np.dot(X, prev))))) + np.dot(pg.transpose(), (
                                    new - prev)) + 0.5 * lam * np.square(linalg.norm(prev - new)):
            return False
        else:
            return True

    def estimatingBeta(self, s11, s22, t):
        '''
        :return: None (betas returned as private variable)
        '''
        iter = 1
        L = np.zeros_like(self.b1)
        while iter < self.maxADMMIter and self.checkConvergence():
            if not self.quiet:
                print '\t\t\tADMM Iteration', iter
            iter += 1
            self.b1 = self.solveBeta(self.X1, self.y1, self.b1, self.b2, L, s22 / (2 * t), 1, self.lam1)
            self.b2 = self.solveBeta(self.X2, self.y2, self.b2, self.b1, -L, s11 / (2 * t), -1, self.lam2)
            L = L + 2 * self.rho * (self.b1 - self.b2)

    def mainCostFunction(self, s11, s22, t):
        tmp = 0
        c1 = s22 / (2 * t)
        tmp += c1 * np.sum(np.square(self.y1 - np.dot(self.X1, self.bpg)).transpose()) * (1.0 / self.X1.shape[0])
        c2 = s11 / (2 * t)
        tmp += c2 * np.sum(np.square(self.y2 - np.dot(self.X2, self.bpg)).transpose()) * (1.0 / self.X2.shape[0])
        tmp += self.lam1 * linalg.norm(self.b1, ord=1) + self.lam2 * linalg.norm(self.b2, ord=1)
        tmp += np.log(t)
        tmp += 2 * self.rho * linalg.norm(self.b1 - self.b2)
        return tmp

    def fit(self, X1, y1, X2, y2, S1=None, U1=None, S2=None, U2=None):
        self.convergenceCheckerCount = 0

        X01 = np.ones(len(y1)).reshape(len(y1), 1)
        X1 = np.hstack([X1, X01])
        X02 = np.ones(len(y2)).reshape(len(y2), 1)
        X2 = np.hstack([X2, X02])

        self.X1 = X1
        self.X2 = X2
        self.y1 = y1
        self.y2 = y2

        [n1, p1] = self.X1.shape
        [n2, p2] = self.X2.shape
        assert p1 == p2

        p = p1
        self.b1 = np.zeros([p, 1])
        self.b2 = np.zeros([p, 1])
        self.b1Prev = None
        self.b2Prev = None

        if not self.quiet:
            print 'Fitting variance'
        tu1, tu2, se1, se2 = self.estimatingVariance(S1, U1, S2, U2)

        iter = 1
        if not self.quiet:
            print 'Running ...'
        while iter < self.maxIter and self.checkConvergence():
            if not self.quiet:
                print '\tIteration:', iter
            iter += 1
            if not self.quiet:
                print '\t\tCalculating Sigma'
            s11, s22, t = self.calculatingSigmaT(tu1, tu2, se1, se2)
            if not self.quiet:
                print '\t\tEstimating Beta'
            self.estimatingBeta(s11, s22, t)

    def getBeta1(self):
        self.b1 = self.b1.reshape(self.b1.shape[0])
        return self.b1[:-1]

    def getBeta2(self):
        self.b2 = self.b2.reshape(self.b2.shape[0])
        return self.b2[:-1]

    def predict1(self, X):
        X0 = np.ones(X.shape[0]).reshape(X.shape[0], 1)
        X = np.hstack([X, X0])
        if not self.logistic:
            return np.dot(X, self.b1)
        else:
            t = 1. / (1 + np.exp(-np.dot(X, self.b1)))
            y = np.zeros_like(t)
            y[t>=0.5] = 1
            return y

    def predict2(self, X):
        X0 = np.ones(X.shape[0]).reshape(X.shape[0], 1)
        X = np.hstack([X, X0])
        if not self.logistic:
            return np.dot(X, self.b2)
        else:
            t = 1. / (1 + np.exp(-np.dot(X, self.b2)))
            y = np.zeros_like(t)
            y[t>=0.5] = 1
            return y
