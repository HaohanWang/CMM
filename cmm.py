# Main file for usage of CMM (Couple Mixed Model)
# If you find this work useful, please cite:
# Wang H., Liu M., Pei F., Vanyukov MM., Bahar I., Wu W. and Xing EP.
# Coupled Mixed Model for Joint GWAS of Pairs of Complex Disorders with Independently Collected Datasets
# Contact: {haohanw,weiwu2,epxing}@cs.cmu.edu

import sys
import numpy as np
from utility.dataLoader import FileReader
from utility.matchSNPList import matchSNPs
from model.CMM import CMM

def printOutHead(out): out.write("\t".join(["RANK", "SNP_ID", "EFFECT_SIZE_ABS"]) + "\n")

def outputResult(out, rank, id, beta):
    out.write("\t".join([str(x) for x in [rank, id, beta]]) + "\n")


from optparse import OptionParser, OptionGroup

usage = """usage: %prog [options] --file1 fileName1 --file2 fileName2
This program provides the basic usage to CMM, e.g:
python cmm.py --file1 data/mice1.plink --file2 data/mice2.plink
	    """
parser = OptionParser(usage=usage)

dataGroup = OptionGroup(parser, "Data Options")
modelGroup = OptionGroup(parser, "Model Options")

## data options
dataGroup.add_option("-t", dest='fileType', default='plink', help="choices of input file type")
dataGroup.add_option("--file1", dest='fileName1', help="name of the input file for one data set")
dataGroup.add_option("--file2", dest='fileName2', help="name of the input file for the other data set")

## model options
modelGroup.add_option("--lambda", dest="lmbd", default=None,
                      help="the weight of the penalizer. If neither lambda or snum is given, cross validation will be run.")
modelGroup.add_option("--rho", dest="rho", default=1e10,
                      help="the weight of the penalizer for common support. If neither lambda or snum is given, cross validation will be run.")
modelGroup.add_option("--snum", dest="snum", default=None,
                      help="the number of targeted variables the model selects. If neither lambda or snum is given, cross validation will be run. "
                           "We highly recommend users to specify the numbers to select, "
                           "then the program will run to make sure there are at least some common snps to be selected for these phenotypes.")
modelGroup.add_option('-q', action='store_true', dest='quiet', default=False, help='Run in quiet mode')
modelGroup.add_option('-m', action='store_true', dest='missing', default=False,
                      help='Run without missing genotype imputation')

## advanced options
parser.add_option_group(dataGroup)
parser.add_option_group(modelGroup)

(options, args) = parser.parse_args()

learningRate = 1e-5

min_lambda_default = 1e-5
max_lambda_default = 1e5
minFactor = 0.6
maxFactor = 1.5
maxIteration = 50
patience = 5


def cross_val_score(clf, X1, y1, X2, y2, cv=5, learningRate=1):
    scores = []
    [n1, p] = X1.shape
    [n2, p] = X2.shape
    b1 = n1 / cv
    b2 = n2 / cv
    for i in range(cv):
        ind1 = np.arange(b1) + b1 * i
        Xtr1 = np.delete(X1, ind1, axis=0)
        ytr1 = np.delete(y1, ind1, axis=0)
        Xte1 = X1[ind1, :]
        yte1 = y1[ind1]

        ind2 = np.arange(b2) + b2 * i
        Xtr2 = np.delete(X2, ind2, axis=0)
        ytr2 = np.delete(y2, ind2, axis=0)
        Xte2 = X1[ind2, :]
        yte2 = y1[ind2]


        clf.setLearningRate(learningRate)
        clf.fit(Xtr1, ytr1, Xtr2, ytr2)
        ypr1 = clf.predict1(Xte1)
        ypr2 = clf.predict2(Xte2)
        s = np.mean(np.square(ypr1 - yte1)) + np.mean(np.square(ypr2 - yte2))
        scores.append(s)
    return scores


def crossValidation(clf, X1, y1, X2, y2, learningRate):
    minError = np.inf
    minLam = 0
    minRho = 0
    for r in range(-10, 10):
        rho = 10.0**r
        for i in range(-10, 10):
            lam = 10.0**i
            clf.setLambda(lam)
            clf.setRho(rho)
            scores = cross_val_score(clf, X1, y1, X2, y2, cv=5, learningRate=learningRate)
            score = np.mean(np.abs(scores))
            if score < minError:
                minError = score
                minLam = lam
                minRho = rho
    clf.setLambda(minLam)
    clf.setRho(minRho)
    clf.setLearningRate(learningRate)
    clf.fit(X1, y1, X2, y2)
    beta1 = clf.getBeta1()
    beta2 = clf.getBeta1()
    return beta1, beta2

def binarySearchJoint_DoubleLambda(model, X1, y1, X2, y2,
                                   select_num, lrScale, quiet=True):
    betaM1 = np.zeros([X1.shape[1]])
    betaM2 = np.zeros([X2.shape[1]])


    min_lambda1 = min_lambda_default
    max_lambda1 = max_lambda_default
    min_lambda2 = min_lambda_default
    max_lambda2 = max_lambda_default

    stuckCount = 1
    previousC = -1
    iteration = 0

    while ((min_lambda1 < max_lambda1) or (min_lambda2 < max_lambda2)) and iteration < maxIteration:
        iteration += 1
        breakFlag = False
        if min_lambda1 < max_lambda1:
            lmbd1 = np.exp((np.log(min_lambda1) + np.log(max_lambda1)) / 2.0)
        if min_lambda2 < max_lambda2:
            lmbd2 = np.exp((np.log(min_lambda2) + np.log(max_lambda2)) / 2.0)

        if not quiet:
            print "\t\tIter:{}\tlambda1:{}\tlambda2:{}".format(iteration, lmbd1, lmbd2),
        model.setLambda1(lmbd1)
        model.setLambda2(lmbd2)
        model.setLearningRate(learningRate * lrScale)  # learning rate must be set again every time we run it.
        model.fit(X1, y1, X2, y2)
        beta1 = model.getBeta1()
        beta2 = model.getBeta2()

        c1 = len(np.where(np.abs(beta1) > 0)[0])
        c2 = len(np.where(np.abs(beta2) > 0)[0])
        c0 = len(np.where(np.minimum(np.abs(beta1), np.abs(beta2))!=0)[0])

        if not quiet:
            print "# Chosen:{} and {} with common {}".format(c1, c2, c0)

        if c1 > select_num * minFactor and c2 > select_num * minFactor and c0 == 0:
            return beta1, beta2

        if c1 < select_num * minFactor:  # Regularizer too strong
            max_lambda1 = lmbd1
        elif c1 > select_num * maxFactor:  # Regularizer too weak
            min_lambda1 = lmbd1
        else:
            breakFlag = True
        if c2 < select_num * minFactor:  # Regularizer too strong
            max_lambda2 = lmbd2
        elif c2 > select_num * maxFactor:  # Regularizer too weak
            min_lambda2 = lmbd2
        else:
            if breakFlag:
                betaM1 = beta1
                betaM2 = beta2
                break
        betaM1 = beta1
        betaM2 = beta2
        if c1 + c2 == previousC:
            stuckCount += 1
        else:
            previousC = c1 + c2
            stuckCount = 1
        if stuckCount > patience:
            if not quiet:
                print 'Run out of patience'
            break

    return betaM1, betaM2

def binarySearch_Rho(model, X1, y1, X2, y2, snum, quiet):
    minRho = 1e0
    maxRho = 1e15
    maxMainIteration = 10
    mainIteration = 0

    Beta1 = None
    Beta2 = None

    while mainIteration < maxMainIteration and minRho < maxRho:
        rho = np.exp((np.log(minRho) + np.log(maxRho)) / 2.0)
        if not quiet:
            print '-------------------------------'
            print 'RHO:', rho

        model.setRho(rho)
        beta1, beta2 = binarySearchJoint_DoubleLambda(model, X1=X1, y1=y1, X2=X2, y2=y2, select_num=snum, lrScale=1e-3, quiet=quiet)

        ind1 = np.where(beta1!=0)[0].tolist()
        ind2 = np.where(beta2!=0)[0].tolist()

        if len(ind1) + len(ind2) > snum*maxFactor*100:
            maxRho = rho
            continue
        elif len(ind1) == 0 and len(ind2) == 0:
            maxRho = rho
            continue
        else:
            m = []
            for a in ind1:
                if a in ind2:
                    m.append(a)
            if not quiet:
                print 'SNPs in common', len(m)

            if len(m) < 1:
                minRho = rho
                continue
            else:
                Beta1 = beta1
                Beta2 = beta2
                break
    return Beta1, Beta2

fileType = 0
IN = None

if len(args) != 0:
    parser.print_help()
    sys.exit()

outFile1 = options.fileName1 + '.output'
outFile2 = options.fileName2 + '.output'

print 'Running ... '

reader = FileReader(fileName=options.fileName1, fileType=options.fileType, imputation=(not options.missing))
X1, Y1, Xname1 = reader.readFiles()

reader = FileReader(fileName=options.fileName2, fileType=options.fileType, imputation=(not options.missing))
X2, Y2, Xname2 = reader.readFiles()

X1, Xname1, X2, Xname2 = matchSNPs(X1, Xname1, X2, Xname2)

model = CMM(quiet=options.quiet)

print 'Computation starts ... '

if options.snum is not None:  # select for a fix number of variable
    snum = float(options.snum)
    Beta1, Beta2 = binarySearch_Rho(model, X1, Y1, X2, Y2, snum, quiet=options.quiet)
elif options.lmbd is not None:
    lmbd = float(options.lmbd)
    model.setLambda(lmbd)
    model.setRho(float(options.rho))
    model.setLearningRate(learningRate)  # learning rate must be set again every time we run it.
    model.fit(X1, Y1, X2, Y2)
    Beta1 = model.getBeta1()
    Beta2 = model.getBeta2()
else:
    Beta1, Beta2 = crossValidation(model, X1, Y1, X2, Y2, learningRate)

ind1 = np.where(Beta1 != 0)[0]
bs1 = Beta1[ind1].tolist()
xname1 = []
for i in ind1:
    xname1.append(i)

beta_name1 = zip(bs1, xname1)
bn = sorted(beta_name1)
bn.reverse()

out1 = open(outFile1, 'w')
printOutHead(out1)

for i in range(len(bn)):
    outputResult(out1, i + 1, bn[i][1], bn[i][0])

out1.close()

ind2 = np.where(Beta2 != 0)[0]
bs2 = Beta2[ind2].tolist()
xname2 = []
for i in ind2:
    xname1.append(i)

beta_name2 = zip(bs2, xname2)
bn = sorted(beta_name2)
bn.reverse()

out2 = open(outFile2, 'w')
printOutHead(out2)

for i in range(len(bn)):
    outputResult(out2, i + 1, bn[i][1], bn[i][0])

out1.close()

print '\nComputation ends normally, check the output file at ', outFile1, outFile2
