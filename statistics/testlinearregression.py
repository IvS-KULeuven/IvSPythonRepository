"""
Unit tests covering linearregression.py

Author: Joris De Ridder
"""


import unittest
from math import pi, sqrt
import numpy as np
from .linearregression import LinearModel, PolynomialModel, HarmonicModel, LinearFit



class LinearModelTestCase(unittest.TestCase):

    """
    Test the LinearModel class and its methods.
    """

    def setUp(self):

       # Generate a simple noiseless dataset

       self.nObservations = 10
       self.nParameters = 2
       time = np.arange(self.nObservations, dtype=np.double)
       self.observations = 3.0 + 2.0 * time**2
       self.regressorList = [np.ones_like(time), time**2]
       self.regressorNames = ["1", "t^2"]
       self.designMatrix = np.column_stack(self.regressorList)
       self.RegressorList2 = [time**3]
       self.RegressorNames2 = ["t^3"]
       self.newnParameters = 3
       self.newModelMatrix = np.vstack(self.regressorList + self.RegressorList2).T


    def tearDown(self):
       pass


    def testInitLinearModel(self):

       self.assertRaises(ValueError, lambda : LinearModel(self.regressorList, ["only 1 name"]))

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertTrue(linearModel.nParameters() == self.nParameters)
       self.assertTrue(linearModel.regressorNames() == self.regressorNames)
       self.assertTrue(linearModel.nObservations() == self.nObservations)


    def testContainsRegressorName(self):

        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertTrue(linearModel.containsRegressorName("t^2"))
        self.assertFalse(linearModel.containsRegressorName("x"))


    def testSomeRegressorNameContains(self):

        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertTrue(linearModel.someRegressorNameContains("^2"))
        self.assertFalse(linearModel.someRegressorNameContains("**2"))


    def testWithIntercept(self):

        # Test a model with an intercept

        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertTrue(linearModel.withIntercept())

        # Test a model without intercept

        linearModel = LinearModel(self.RegressorList2, self.RegressorNames2)
        self.assertFalse(linearModel.withIntercept())


    def testModelMatrix(self):

       # With a list of regressor arrays as input

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertTrue(isinstance(linearModel.designMatrix(), np.ndarray))
       self.assertTrue(linearModel.designMatrix().dtype == np.double)
       self.assertTrue(linearModel.designMatrix().shape == (self.nObservations, self.nParameters))
       self.assertTrue(np.alltrue(linearModel.designMatrix() == self.designMatrix))

       # With a design matrix as input

       linearModel = LinearModel(self.designMatrix, self.regressorNames)
       self.assertTrue(isinstance(linearModel.designMatrix(), np.ndarray))
       self.assertTrue(linearModel.designMatrix().dtype == np.double)
       self.assertTrue(linearModel.designMatrix().shape == (self.nObservations, self.nParameters))
       self.assertTrue(np.alltrue(linearModel.designMatrix() == self.designMatrix))


    def testPseudoInverse(self):

        expectedPseudoInverse = np.array([[ 0.21264822,  0.20869565,  0.19683794,  0.1770751,   0.14940711,
                                            0.11383399,  0.07035573,  0.01897233, -0.04031621, -0.10750988],
                                          [-0.00395257, -0.00381388, -0.00339782, -0.00270439, -0.00173358,
                                           -0.0004854,   0.00104015,  0.00284308,  0.00492338,  0.00728105]])
        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertTrue(np.allclose(linearModel.pseudoInverse(), expectedPseudoInverse,
                                    rtol=1.0e-6, atol=1.e-08))


    def testStandardCovarianceMatrix(self):

       expectedStandardCovarianceMatrix = np.array([[ 2.12648221e-01, -3.95256917e-03],
                                                    [-3.95256917e-03,  1.38686638e-04]])
       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertTrue(np.allclose(linearModel.standardCovarianceMatrix(),
                                   expectedStandardCovarianceMatrix, rtol=1.0e-6, atol=1.e-08))


    def testConditionNumber(self):

        expectedConditionNumber = 57.120830024
        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertAlmostEqual(linearModel.conditionNumber(), expectedConditionNumber, places = 5)


    def testSingularValues(self):

        expectedSingularValues = np.array([123.84788663, 2.16817379])
        linearModel = LinearModel(self.regressorList, self.regressorNames)
        self.assertTrue(np.allclose(linearModel.singularValues(), expectedSingularValues,
                                    rtol=1.0e-6, atol=1.e-08))


    def testHatMatrix(self):

       expectedHatMatrix = np.array([[ 0.21264822,  0.20869565,  0.19683794,  0.1770751 ,  0.14940711,
                                       0.11383399,  0.07035573,  0.01897233, -0.04031621, -0.10750988],
                                     [ 0.20869565,  0.20488177,  0.19344012,  0.17437071,  0.14767353,
                                       0.11334859,  0.07139588,  0.02181541, -0.03539283, -0.10022883],
                                     [ 0.19683794,  0.19344012,  0.18324665,  0.16625754,  0.14247278,
                                       0.11189238,  0.07451633,  0.03034464, -0.0206227 , -0.07838569],
                                     [ 0.1770751 ,  0.17437071,  0.16625754,  0.15273559,  0.13380487,
                                       0.10946536,  0.07971708,  0.04456002,  0.00399418, -0.04198045],
                                     [ 0.14940711,  0.14767353,  0.14247278,  0.13380487,  0.12166979,
                                       0.10606754,  0.08699813,  0.06446155,  0.0384578 ,  0.00898689],
                                     [ 0.11383399,  0.11334859,  0.11189238,  0.10946536,  0.10606754,
                                       0.10169891,  0.09635948,  0.09004923,  0.08276819,  0.07451633],
                                     [ 0.07035573,  0.07139588,  0.07451633,  0.07971708,  0.08699813,
                                       0.09635948,  0.10780112,  0.12132307,  0.13692532,  0.15460786],
                                     [ 0.01897233,  0.02181541,  0.03034464,  0.04456002,  0.06446155,
                                       0.09004923,  0.12132307,  0.15828306,  0.2009292 ,  0.24926149],
                                     [-0.04031621, -0.03539283, -0.0206227 ,  0.00399418,  0.0384578 ,
                                       0.08276819,  0.13692532,  0.2009292 ,  0.27477983,  0.35847722],
                                     [-0.10750988, -0.10022883, -0.07838569, -0.04198045,  0.00898689,
                                       0.07451633,  0.15460786,  0.24926149,  0.35847722,  0.48225504]])

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertTrue(np.allclose(linearModel.hatMatrix(), expectedHatMatrix, rtol=1.0e-6, atol=1.e-08))


    def testTraceHatMatrix(self):

        linearModel = LinearModel(self.regressorList, self.regressorNames)
        expectedTrace = 2.0
        self.assertAlmostEqual(linearModel.traceHatMatrix(), expectedTrace, places=7)


    def testDegreesOfFreedom(self):

       expectedDegreesOfFreedom = 8
       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertEqual(linearModel.degreesOfFreedom(), expectedDegreesOfFreedom)


    def testAddLinearModel(self):

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       linearModel2 = LinearModel(self.RegressorList2, self.RegressorNames2)
       linearModel = linearModel + linearModel2
       self.assertTrue(linearModel.nParameters() == self.newnParameters)
       self.assertTrue(np.alltrue(linearModel.designMatrix() == self.newModelMatrix))
       self.assertTrue(linearModel.regressorNames() == self.regressorNames + self.RegressorNames2)
       self.assertTrue(linearModel._covMatrixObserv == None)

       self.assertRaises(TypeError, lambda lm : lm + 1.0, linearModel)
       self.assertRaises(ValueError, lambda x,y :  LinearModel(x,y) + LinearModel(x[:-1,:],y),
                         self.designMatrix, self.regressorNames)


    def testCopyLinearModel(self):

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       linearModel2 = linearModel.copy()
       self.assertTrue(linearModel2.nParameters() == self.nParameters)
       self.assertTrue(np.alltrue(linearModel2.designMatrix() == self.designMatrix))
       self.assertTrue(linearModel2.regressorNames() == self.regressorNames)
       self.assertTrue(id(linearModel2.designMatrix()) != id(linearModel.designMatrix()))
       self.assertTrue(id(linearModel2.regressorNames()) != id(linearModel.regressorNames()))


    def testFitData(self):

       linearModel = LinearModel(self.regressorList, self.regressorNames)
       self.assertTrue(isinstance(linearModel.fitData(self.observations), LinearFit))





class WeightedLinearModelTestCase(unittest.TestCase):

    """
    Test the LinearModel class with a specified covariance matrix for the
    observations. Test only those method affected by weights.
    """

    def setUp(self):

        # Generate simple noiseless dataset

        self.nObservations = 10
        self.nParameters = 2
        self.time = np.arange(self.nObservations, dtype=np.double)
        self.regressorList = [np.ones_like(self.time), self.time]
        self.unweightedDesignMatrix = np.column_stack(self.regressorList)
        self.regressorNames = ["1", "t"]

        # First covariance matrix: different weights, no correlations

        self.covMatrixObserv1 = np.diag(0.1 + 0.1 * np.arange(self.nObservations))

        # Second covariance matrix: different weights and correlations
        # Correlation coefficient: 0.6

        self.covMatrixObserv2 = np.diag(0.1 + 0.1 * np.arange(self.nObservations))
        for i in range(self.nObservations):
            for j in range(self.nObservations):
                if i>=j: continue
                self.covMatrixObserv2[i,j] = 0.6 * sqrt(self.covMatrixObserv2[i,i] * self.covMatrixObserv2[j,j])
                self.covMatrixObserv2[j,i] = self.covMatrixObserv2[i,j]


    def tearDown(self):
       pass


    def testInitWeightedLinearModel(self):

        # The covariance matrix of the observations should be 1) a numpy array
        # 2) with the proper dimensions, and 3) with the correct size.

        self.assertRaises(TypeError, lambda x,y,z : LinearModel(x,y,z),             \
                          self.regressorList, self.regressorNames, [1,2,3])
        self.assertRaises(TypeError, lambda x,y,z : LinearModel(x,y,z),             \
                          self.regressorList, self.regressorNames, np.arange(10))
        self.assertRaises(TypeError, lambda x,y,z : LinearModel(x,y,z),             \
                          self.regressorList, self.regressorNames,                  \
                          np.arange(self.nObservations*5).reshape(self.nObservations,5))

        linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv1, regressorsAreWeighted=True)
        self.assertTrue(isinstance(linearModel._covMatrixObserv, np.ndarray))
        self.assertTrue(linearModel._covMatrixObserv.dtype == np.double)
        self.assertTrue(linearModel._covMatrixObserv.shape == (self.nObservations, self.nObservations))
        self.assertTrue(np.alltrue(linearModel._covMatrixObserv == self.covMatrixObserv1))

        linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv1, regressorsAreWeighted=False)
        self.assertTrue(isinstance(linearModel._covMatrixObserv, np.ndarray))
        self.assertTrue(linearModel._covMatrixObserv.dtype == np.double)
        self.assertTrue(linearModel._covMatrixObserv.shape == (self.nObservations, self.nObservations))
        self.assertTrue(np.alltrue(linearModel._covMatrixObserv == self.covMatrixObserv1))


    def testWeightedModelMatrix(self):

       linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2, regressorsAreWeighted=True)
       self.assertTrue(np.alltrue(linearModel.designMatrix() == self.unweightedDesignMatrix))

       linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2, regressorsAreWeighted=False)
       self.assertFalse(np.alltrue(linearModel.designMatrix() == self.unweightedDesignMatrix))

       expectedWeightedDesignMatrix = np.array([[ 3.16227766e+00, 4.73643073e-15],
                                                [ 4.23376727e-01, 2.79508497e+00],
                                                [-2.67843095e-01, 3.79299210e+00],
                                                [-5.45288776e-01, 4.39760881e+00],
                                                [-6.78144517e-01, 4.84809047e+00],
                                                [-7.46997591e-01, 5.21965041e+00],
                                                [-7.83412990e-01, 5.54374745e+00],
                                                [-8.01896636e-01, 5.83605564e+00],
                                                [-8.09868157e-01, 6.10538773e+00],
                                                [-8.11425366e-01, 6.35716086e+00]])

       self.assertTrue(np.allclose(linearModel.designMatrix(), expectedWeightedDesignMatrix, rtol=1.0e-6, atol=1.e-08))


    def testWeightedHatMatrix(self):

        expectedWeightedHatMatrix = np.array([[ 0.93746979,  0.22625449,  0.05730702, -0.00315076, -0.02629974,
                                               -0.03331944, -0.03243355, -0.02737749, -0.02003319, -0.01142021],
                                              [ 0.22625449,  0.08788299,  0.05898899,  0.05159613,  0.05137251,
                                                0.05410201,  0.05817441,  0.0628748 ,  0.0678539 ,  0.07293014],
                                              [ 0.05730702,  0.05898899,  0.06478378,  0.07085639,  0.07671941,
                                                0.08229333,  0.08758368,  0.09261537,  0.09741573,  0.10200995],
                                              [-0.00315076,  0.05159613,  0.07085639,  0.08238505,  0.09090112,
                                                0.09788464,  0.10395253,  0.10941096,  0.11443131,  0.11911848],
                                              [-0.02629974,  0.05137251,  0.07671941,  0.09090112,  0.1008532 ,
                                                0.10872302,  0.11539093,  0.12128539,  0.12664119,  0.13159879],
                                              [-0.03331944,  0.05410201,  0.08229333,  0.09788464,  0.10872302,
                                                0.11723346,  0.12440767,  0.13072689,  0.13645396,  0.14174555],
                                              [-0.03243355,  0.05817441,  0.08758368,  0.10395253,  0.11539093,
                                                0.12440767,  0.13203013,  0.13875766,  0.14486348,  0.15051078],
                                              [-0.02737749,  0.0628748 ,  0.09261537,  0.10941096,  0.12128539,
                                                0.13072689,  0.13875766,  0.14587641,  0.15235719,  0.15836442],
                                              [-0.02003319,  0.0678539 ,  0.09741573,  0.11443131,  0.12664119,
                                                0.13645396,  0.14486348,  0.15235719,  0.15920448,  0.16556802],
                                              [-0.01142021,  0.07293014,  0.10200995,  0.11911848,  0.13159879,
                                                0.14174555,  0.15051078,  0.15836442,  0.16556802,  0.17228071]])

        linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2)
        self.assertTrue(np.allclose(linearModel.hatMatrix(), expectedWeightedHatMatrix, rtol=1.0e-6, atol=1.e-08))


    def testWeightedDegreesOfFreedom(self):

        expectedDegreesOfFreedom = 8
        linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2)
        self.assertAlmostEqual(linearModel.degreesOfFreedom(), expectedDegreesOfFreedom, places=6)


    def testAddWeigthedLinearModel(self):

        regressorList2 = [self.time**2]
        regressorNames2 = ["t^2"]

        self.assertRaises(ValueError, lambda v,w,x,y,z :  LinearModel(v,w) + LinearModel(x,y,z),
                          self.regressorList, self.regressorNames,
                          regressorList2, regressorNames2, self.covMatrixObserv1)

        self.assertRaises(ValueError, lambda u,v,w,x,y,z :  LinearModel(u,v,w) + LinearModel(x,y,z),
                          self.regressorList, self.regressorNames, self.covMatrixObserv1,
                          regressorList2, regressorNames2, self.covMatrixObserv2)

        linearModel1 = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2)
        linearModel2 = LinearModel(regressorList2, regressorNames2, self.covMatrixObserv2)
        linearModel = linearModel1 + linearModel2
        self.assertTrue(np.alltrue(linearModel._covMatrixObserv == self.covMatrixObserv2))


    def testCopyWeightedLinearModel(self):

       linearModel = LinearModel(self.regressorList, self.regressorNames, self.covMatrixObserv2, regressorsAreWeighted=False)
       linearModel2 = linearModel.copy()
       self.assertTrue(np.alltrue(linearModel2._covMatrixObserv == self.covMatrixObserv2))
       self.assertTrue(id(linearModel2._covMatrixObserv) != id(linearModel._covMatrixObserv))





class PolynomialModelTestCase(unittest.TestCase):

    """
    Test the PolynomialModel class and its methods.
    """


    def setUp(self):

       # Generate a simple noiseless dataset

       self.nObservations = 10
       self.time = np.arange(self.nObservations, dtype=np.double)


    def tearDown(self):
       pass


    def testInitPolynomialModel(self):

       # Define the input, and derive what the output should be

       orderList = [0,2]
       nParameters = len(orderList)
       covariateName = "t"
       regressorNames = ["1", "t^2"]
       designMatrix = np.empty((self.nObservations, nParameters))
       for n in range(len(orderList)):
           designMatrix[:,n] = self.time**orderList[n]

       # Assert if the output is what it should be

       polyModel = PolynomialModel(self.time, covariateName, orderList)
       self.assertTrue(polyModel.nParameters() == nParameters)
       self.assertTrue(polyModel.regressorNames() == regressorNames)
       self.assertTrue(polyModel.nObservations() == self.nObservations)
       self.assertTrue(np.alltrue(polyModel.designMatrix() == designMatrix))





class WeightedPolynomialModelTestCase(unittest.TestCase):

    """
    Test if a covariance matrix can be specified for a polynomial model
    """


    def setUp(self):

       # Generate a simple noiseless dataset

       self.nObservations = 10
       self.time = np.arange(self.nObservations, dtype=np.double)


    def tearDown(self):
       pass


    def testInitPolynomialModel(self):

       # Defining a polynomial model with no covariance matrix should be the
       # same as one with a diagonal covariance matrix with ones on the diagonal.

       orderList = [0,2]
       nParameters = len(orderList)
       covariateName = "t"
       regressorNames = ["1", "t^2"]

       polyModel1 = PolynomialModel(self.time, covariateName, orderList)
       polyModel2 = PolynomialModel(self.time, covariateName, orderList, np.diag(np.ones(self.nObservations)))

       self.assertTrue(np.alltrue(polyModel1.designMatrix() == polyModel2.designMatrix()))





class HarmonicModelTestCase(unittest.TestCase):

    """
    Test the HarmonicModel class and its methods.
    """


    def setUp(self):

       # Generate a simple noiseless dataset

       self.nObservations = 10
       self.time = np.arange(self.nObservations, dtype=np.double)


    def tearDown(self):
       pass


    def testInitHarmonicModel(self):

       # Define the input, and derive what the output should be

       freqList = [2.1, 3.4]
       freqNameList = ["f_1", "f_2"]
       maxNharmonics = 2
       nParameters = 2*maxNharmonics*len(freqList)
       covariateName = "t"
       regressorNames = ["sin(2*pi*1*f_1*t)", "cos(2*pi*1*f_1*t)",       \
                         "sin(2*pi*2*f_1*t)", "cos(2*pi*2*f_1*t)",       \
                         "sin(2*pi*1*f_2*t)", "cos(2*pi*1*f_2*t)",       \
                         "sin(2*pi*2*f_2*t)", "cos(2*pi*2*f_2*t)"]
       designMatrix = np.empty((self.nObservations, nParameters))
       runningIndex = 0
       for n in range(len(freqList)):
            for k in range(1,maxNharmonics+1):
               designMatrix[:,runningIndex] = np.sin(2*pi*k*freqList[n]*self.time)
               designMatrix[:,runningIndex+1] = np.cos(2*pi*k*freqList[n]*self.time)
               runningIndex += 2

       # Assert if the output is what it should be

       harmonicModel = HarmonicModel(self.time, covariateName, freqList, freqNameList,maxNharmonics)
       self.assertTrue(harmonicModel.nParameters() == nParameters)
       self.assertTrue(harmonicModel.regressorNames() == regressorNames)
       self.assertTrue(harmonicModel.nObservations() == self.nObservations)
       self.assertTrue(np.alltrue(harmonicModel.designMatrix() == designMatrix))





class WeightedHarmonicModelTestCase(unittest.TestCase):

    """
    Test if a covariance matrix can be specified for a Harmonic model
    """


    def setUp(self):

       # Generate a simple noiseless dataset

       self.nObservations = 10
       self.time = np.arange(self.nObservations, dtype=np.double)


    def tearDown(self):
       pass


    def testInitHarmonicModel(self):

       # Defining a harmonic model with no covariance matrix should be the
       # same as one with a diagonal covariance matrix with ones on the diagonal.

       freqList = [2.1]
       freqNameList = ["f_1"]
       maxNharmonics = 1
       nParameters = 2*maxNharmonics*len(freqList)
       covariateName = "t"
       regressorNames = ["sin(2*pi*1*f_1*t)", "cos(2*pi*1*f_1*t)"]
       covMatrixObserv = np.diag(np.ones(self.nObservations))
       harmonicModel1 = HarmonicModel(self.time, covariateName, freqList, freqNameList, maxNharmonics)
       harmonicModel2 = HarmonicModel(self.time, covariateName, freqList, freqNameList, maxNharmonics, covMatrixObserv)

       self.assertTrue(np.alltrue(harmonicModel1.designMatrix() == harmonicModel2.designMatrix()))





class LinearSubmodelGeneratorTestCase(unittest.TestCase):

    """
    Test the LinearModel.submodels() method
    """

    def setUp(self):

       self.nObservations = 10
       self.x = np.arange(self.nObservations)
       x = self.x   # local shorter alias
       ones= np.ones_like(x)
       self.linearModel1 = LinearModel([ones,x,x**2], ["1", "x", "x^2"])
       self.linearModel2 = LinearModel([ones,x,np.sin(x),np.cos(x)], ["1","x","sin(x)","cos(x)"])


    def tearDown(self):
       pass


    def testSubmodel(self):

        generator = self.linearModel1.submodels(Nmin=3, nested=True, ranks=[1,0,3])
        self.assertTrue(generator.__class__.__name__ == "generator")
        submodel = next(generator)
        self.assertTrue(isinstance(submodel, LinearModel))
        designMatrix = np.column_stack([self.x,np.ones(len(self.x)),self.x**2])
        self.assertTrue(np.alltrue(submodel.designMatrix() == designMatrix))
        self.assertTrue(submodel.regressorNames() == ["x","1","x^2"])


    def testSubmodelGenerator(self):

       names = []
       for model in self.linearModel1.submodels(Nmin=1, nested=False, ranks=None):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1"], ["x"], ["x^2"], ["1","x"], ["1","x^2"], ["x","x^2"], ["1","x","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=1, Nmax=2, nested=False, ranks=None):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1"], ["x"], ["x^2"], ["1","x"], ["1","x^2"], ["x","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=2, nested=False, ranks=None):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1","x"], ["1","x^2"], ["x","x^2"], ["1","x","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=3, nested=True, ranks=None):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1","x","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=3, nested=True, ranks=[1,0,3]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["x","1","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=2, nested=False, ranks=[1,0,3]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["x","1"], ["x","x^2"], ["1","x^2"], ["x","1","x^2"]])

       names = []
       for model in self.linearModel1.submodels(Nmin=2, nested=True, ranks=[1,0,3]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["x","1"], ["x","1","x^2"]])


       names = []
       for model in self.linearModel2.submodels(Nmin=3, nested=True, ranks=[0,1,2,2]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1","x","sin(x)","cos(x)"]])

       names = []
       for model in self.linearModel2.submodels(Nmin=2, nested=True, ranks=[0,1,2,2]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1","x"], ["1","x","sin(x)","cos(x)"]])

       names = []
       for model in self.linearModel2.submodels(Nmin=2, nested=False, ranks=[0,1,2,2]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1","x"], ["1","sin(x)","cos(x)"], ["x","sin(x)","cos(x)"], ["1","x","sin(x)","cos(x)"]])

       names = []
       for model in self.linearModel2.submodels(Nmin=1, Nmax=2, nested=False, ranks=[0,1,2,2]):
          names += [model.regressorNames()]
       self.assertTrue(names == [["1"], ["x"], ["sin(x)","cos(x)"], ["1","x"], ["1","sin(x)","cos(x)"], ["x","sin(x)","cos(x)"]])





class LinearFitTestCase(unittest.TestCase):

    """
    Test the LinearFit class and its methods

    """

    def setUp(self):

       # Generate a simple dataset.
       # Noise-values are drawn from a N(0,1), but are truncated.

       self.nObservations = 10
       self.nParameters = 2
       self.time = np.arange(self.nObservations, dtype=np.double)
       self.noise = np.array([-.72, -1.02, .52, .24, -.17, -1.78, .56, .33, .19, 1.03])
       self.observations = 3.0 + 2.0 * self.time**2 + self.noise
       regressorList = [np.ones_like(self.time), self.time**2]
       regressorNames = ["1", "t^2"]
       self.linearModel = LinearModel(regressorList, regressorNames)
       self.linearFit = self.linearModel.fitData(self.observations)


    def tearDown(self):
        pass


    def testLinearFit(self):

        self.assertTrue(isinstance(self.linearFit, LinearFit))


    def testObservations(self):

        # Since the linear model did not have a covariance matrix specified, the
        # decorrelated observations, should be simply the original observations

        self.assertTrue(isinstance(self.linearFit.observations(weighted=True), np.ndarray))
        self.assertEqual(self.linearFit.observations(weighted=True).shape, (self.nObservations,))
        self.assertTrue(np.alltrue(self.observations == self.linearFit.observations(weighted=True)))
        self.assertTrue(np.alltrue(self.observations == self.linearFit.observations(weighted=False)))



    def testRegressionCoefficients(self):

       expectedCoeff = np.array([2.47811858, 2.01543444])
       self.assertTrue(isinstance(self.linearFit.regressionCoefficients(), np.ndarray))
       self.assertEqual(self.linearFit.regressionCoefficients().shape, (self.nParameters,))
       self.assertTrue(np.allclose(self.linearFit.regressionCoefficients(), expectedCoeff, rtol=1.0e-6, atol=1.e-08))


    def testResiduals(self):

       expectedResiduals = np.array([-0.19811858, -0.51355301, 0.98014368, 0.6229715,
                                      0.10493045, -1.64397947, 0.52624173, 0.09559406,
                                     -0.27592247,  0.30169212])
       self.assertTrue(isinstance(self.linearFit.residuals(weighted=False), np.ndarray))
       self.assertEqual(self.linearFit.residuals(weighted=False).shape, (self.nObservations,))
       self.assertTrue(np.allclose(self.linearFit.residuals(weighted=False), expectedResiduals, rtol=1.0e-6, atol=1.e-08))


    def testPredictions(self):

       expectedPredictions = np.array([2.47811858, 4.49355301, 10.53985632, 20.6170285,
                                       34.72506955, 52.86397947, 75.03375827,
                                       101.23440594, 131.46592247, 165.72830788])
       self.assertTrue(isinstance(self.linearFit.predictions(weighted=False), np.ndarray))
       self.assertEqual(self.linearFit.predictions(weighted=False).shape,  (self.nObservations,))
       self.assertTrue(np.allclose(self.linearFit.predictions(weighted=False), expectedPredictions, rtol=1.0e-6, atol=1.e-08))


    def testCovarianceMatrix(self):

       expectedCovMatrix = np.array([[1.28084978e-01, -2.38076167e-03],
                                     [-2.38076167e-03, 8.35354974e-05]])
       self.assertTrue(isinstance(self.linearFit.covarianceMatrix(), np.ndarray))
       self.assertEqual(self.linearFit.covarianceMatrix().shape,  (self.nParameters,self.nParameters))
       self.assertTrue(np.allclose(self.linearFit.covarianceMatrix(), expectedCovMatrix, rtol=1.0e-6, atol=1.e-08))


    def testErrorBars(self):

       expectedErrors = np.array([0.35788962, 0.00913978])
       self.assertTrue(isinstance(self.linearFit.errorBars(), np.ndarray))
       self.assertEqual(self.linearFit.errorBars().shape, (self.nParameters,))
       self.assertTrue(np.allclose(self.linearFit.errorBars(), expectedErrors, rtol=1.0e-6, atol=1.e-08))


    def testConfidenceIntervals(self):

        alpha = 0.05
        expectedLower = np.array([1.65282364, 1.99435808])
        expectedUpper = np.array([3.30341351, 2.0365108 ])
        self.assertEqual(len(self.linearFit.confidenceIntervals(0.05)), 2)
        lower, upper = self.linearFit.confidenceIntervals(0.05)
        self.assertTrue(isinstance(lower, np.ndarray))
        self.assertTrue(isinstance(upper, np.ndarray))
        self.assertEqual(len(lower), self.nParameters)
        self.assertEqual(len(upper), self.nParameters)
        self.assertTrue(np.allclose(lower, expectedLower, rtol=1.0e-6, atol=1.e-08))
        self.assertTrue(np.allclose(upper, expectedUpper, rtol=1.0e-6, atol=1.e-08))


    def testT_values(self):

        expectedTvalues = np.array([6.92425385, 220.51235807])
        self.assertTrue(isinstance(self.linearFit.t_values(), np.ndarray))
        self.assertEqual(self.linearFit.t_values().shape, (self.nParameters,))
        self.assertTrue(np.allclose(self.linearFit.t_values(), expectedTvalues, rtol=1.0e-6, atol=1.e-08))


    def testSumOfSquaredResiduals(self):

        expectedSumSqResiduals = 4.81866162957
        self.assertAlmostEqual(self.linearFit.sumSqResiduals(weighted=False), expectedSumSqResiduals, places=7)


    def testResidualVariance(self):

        expectedResidualVariance = 0.60233270369599989
        self.assertAlmostEqual(self.linearFit.residualVariance(weighted=False), expectedResidualVariance, places=7)


    def testCorrelationMatrix(self):

       expectedCorrelationMatrix = np.array([[1., -0.72783225], [-0.72783225, 1.]])
       self.assertTrue(np.allclose(self.linearFit.correlationMatrix(), expectedCorrelationMatrix, rtol=1.0e-6, atol=1.e-08))


    def testCoefficientOfDetermination(self):

        # Test using a model with an intercept

        expectedCoefficientOfDetermination = 0.99983550516908659
        self.assertAlmostEqual(self.linearFit.coefficientOfDetermination(weighted=False), expectedCoefficientOfDetermination, places=7)

        # Test using a model without an intercept

        linearModel = LinearModel([self.time**2], ["t^2"])
        linearFit = linearModel.fitData(self.observations)
        expectedCoefficientOfDetermination = 0.99948312767720937
        self.assertAlmostEqual(linearFit.coefficientOfDetermination(weighted=False), expectedCoefficientOfDetermination, places=7)


    def testBICvalue(self):

        expectedBICvalue = -0.39313345804942657
        self.assertAlmostEqual(self.linearFit.BICvalue(), expectedBICvalue, places=7)


    def testAICvalue(self):

        expectedAICvalue = 2.6991112629684357
        self.assertAlmostEqual(self.linearFit.AICvalue(), expectedAICvalue, places=7)


    def testThypothesisTest(self):

        regressorList = [np.ones_like(self.time), self.time**2, self.time**3]
        regressorNames = ["1", "t^2", "t^3"]
        linearModel = LinearModel(regressorList, regressorNames)
        linearFit = linearModel.fitData(self.observations)

        alpha = 0.05                                     # significance level
        expectedOutput = np.array([True, True, False])
        nullRejected = linearFit.regressorTtest(alpha)

        self.assertTrue(isinstance(nullRejected, np.ndarray))
        self.assertTrue(nullRejected.dtype == np.bool)
        self.assertEqual(len(nullRejected), len(regressorList))
        self.assertTrue(np.all(nullRejected == expectedOutput))


    def testFstatistic(self):

        # Test on a model with more than one regressor

        expectedFstatistic = 48625.747064132738
        self.assertAlmostEqual(self.linearFit.Fstatistic(weighted=False), expectedFstatistic, places=7)

        # Test on a model with only one regressor

        expectedFstatistic = 17403.423925909559
        linearModel = LinearModel([np.square(self.time)], ["t^2"])
        linearFit = linearModel.fitData(self.observations)
        self.assertAlmostEqual(linearFit.Fstatistic(weighted=False), expectedFstatistic, places=7)


    def testFstatisticTest(self):

        # with tested model the same as the true model

        alpha = 0.05
        expectedOutput = True         # null hypothesis rejected
        self.assertTrue(self.linearFit.FstatisticTest(alpha) == expectedOutput)


        # with a fit on pure noise noise

        alpha = 0.05
        expectedOutput = False        # null hypothesis could not be rejected
        linearFit = self.linearModel.fitData(self.noise)
        self.assertTrue(linearFit.FstatisticTest(alpha) == expectedOutput)



    def testEvaluate(self):

        expectedOutput = np.array([ 52.86397947,  42.8147219 ,  33.8820485 ,  26.06595927,
                                    19.36645422,  13.78353335,   9.31719665,   5.96744412,
                                    3.73427577,   2.6176916 ,   2.6176916 ,   3.73427577,
                                    5.96744412,   9.31719665,  13.78353335,  19.36645422,
                                    26.06595927,  33.8820485 ,  42.8147219 ,  52.86397947])

        xnew = np.linspace(-5.0, +5.0, 20)
        newRegressorList = [np.ones_like(xnew), xnew**2]
        output = self.linearFit.evaluate(newRegressorList)
        self.assertTrue(np.allclose(output, expectedOutput, rtol=1.0e-6, atol=1.e-08))





class WeightedLinearFitTestCase(unittest.TestCase):

    """
    Test the LinearFit class with a specified covariance matrix for the
    observations. Test only those method affected by weights.

    """

    def setUp(self):

        # Generate simple noiseless dataset

        self.nObservations = 10
        self.nParameters = 2
        time = np.arange(self.nObservations, dtype=np.double)
        self.observations1 = 1.0 + 0.3 * time
        self.observations2 = 1.0 + 0.3 * time
        self.regressorList = [np.ones_like(time), time]
        self.regressorNames = ["1", "t"]
        self.designMatrix = np.column_stack(self.regressorList)

        # First covariance matrix: different weights, no correlations

        covMatrixObserv1 = np.diag(0.1 + 0.1 * np.arange(self.nObservations))

        # Second covariance matrix: different weights and correlations
        # Correlation coefficient: 0.6

        covMatrixObserv2 = np.diag(0.1 + 0.1 * np.arange(self.nObservations))
        for i in range(self.nObservations):
            for j in range(self.nObservations):
                if i>=j: continue
                covMatrixObserv2[i,j] = 0.6 * sqrt(covMatrixObserv2[i,i] * covMatrixObserv2[j,j])
                covMatrixObserv2[j,i] = covMatrixObserv2[i,j]

        # Add noise according to the different covariance matrices

        noise1 = np.array([-0.26944792, -0.56542802,  0.62106263, -0.03113657,  0.98090236,      \
                            0.02678669,  1.2237701 , -0.50112787, -0.47742454,  1.16351356])
        noise2 = np.array([ 0.24484502, -0.22979797,  0.40639882,  0.12137103,  0.20694025,      \
                            0.68952746, -0.30300402, -0.11136982,  0.3549814 ,  0.20528704])

        self.observations1 += noise1
        self.observations2 += noise2

        # Instantiate the (weighted) linear models and their linear fits

        self.linearModel1 = LinearModel(self.regressorList, self.regressorNames, covMatrixObserv1)
        self.linearModel2 = LinearModel(self.regressorList, self.regressorNames, covMatrixObserv2)
        self.linearFit1 = self.linearModel1.fitData(self.observations1)
        self.linearFit2 = self.linearModel2.fitData(self.observations2)


    def tearDown(self):
       pass


    def testDecorrelatedObservations(self):

        expectedObservations = np.array([3.9365456,  0.03889641, 1.73885587, 0.65979961, 0.82899163,
                                         1.73402234, -0.1799988, 0.37289876, 1.25537791,  1.05209144])

        self.assertTrue(isinstance(self.linearFit2.observations(weighted=True), np.ndarray))
        self.assertEqual(self.linearFit2.observations(weighted=True).shape, (self.nObservations,))
        self.assertTrue(np.allclose(expectedObservations, self.linearFit2.observations(weighted=True), rtol=1.0e-6, atol=1.e-08))


    def testWeightedRegressionCoefficients(self):

        # A test with a diagonal covariance matrix

        expectedWeightedCoeff = np.array([0.76576446, 0.40030724])
        self.assertTrue(np.allclose(self.linearFit1.regressionCoefficients(), expectedWeightedCoeff, rtol=1.0e-6, atol=1.e-08))

        # A test with a non-diagonal covariance matrix

        expectedWeightedCoeff = np.array([1.1623421, 0.3040605])
        self.assertTrue(np.allclose(self.linearFit2.regressionCoefficients(), expectedWeightedCoeff, rtol=1.0e-6, atol=1.e-08))


    def testResiduals(self):

        # A test with a diagonal covariance matrix

        expectedResiduals = np.array([-0.03521238, -0.43149972, 0.65468369, -0.09782275, 0.81390894, -0.24051397, 0.85616220, -0.96904301, -1.04564692, 0.49498394])
        self.assertTrue(np.allclose(self.linearFit1.residuals(weighted=False), expectedResiduals, rtol=1.0e-6, atol=1.e-08))

        # A test with a non-diagonal covariance matrix

        expectedResiduals = np.array([0.082502896, -0.396200554, 0.235935777, -0.053152473, 0.028356288, 0.506883039, -0.489708901, -0.302135160, 0.160155600, 0.006400781])
        self.assertTrue(np.allclose(self.linearFit2.residuals(weighted=False), expectedResiduals, rtol=1.0e-6, atol=1.e-08))












suite = [unittest.TestLoader().loadTestsFromTestCase(LinearModelTestCase)]
suite += unittest.TestLoader().loadTestsFromTestCase(WeightedLinearModelTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(PolynomialModelTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(WeightedPolynomialModelTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(HarmonicModelTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(WeightedHarmonicModelTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(LinearSubmodelGeneratorTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(LinearFitTestCase)
suite += unittest.TestLoader().loadTestsFromTestCase(WeightedLinearFitTestCase)

allTests = unittest.TestSuite(suite)
unittest.TextTestRunner(verbosity=2).run(allTests)
