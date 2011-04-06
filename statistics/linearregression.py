"""
Easy linear fitting.
Author: Joris De Ridder

A crash course
==============
The following example will illustrate some of the core functionality of the package.
We first make some "fake" observations I{obs} that were measured as a function of I{x}.

>>> from numpy import *
>>> noise = array([0.44, -0.48, 0.26, -2.00, -0.93, 2.21, -0.57, -2.04, -1.09, 1.53])
>>> x = linspace(0, 5, 10)
>>> obs = 2.0 + 3.0 * exp(x) + noise          # our "observations"

The model we want to fit is y = a_0 + a_1 * exp(x) = a_0 * 1 + a_1 * exp(x):

>>> myModel = LinearModel([ones(10), exp(x)], ["1", "exp(x)"])
>>> print myModel
Model: y = a_0 + a_1 * exp(x)
Expected number of observations: 10

Now we fit the observations:

>>> myFit = myModel.fitData(obs)
>>> myFit.summary()
Model: y = a_0 + a_1 * exp(x)
Regression coefficients:
    a_0 = 1.519644 +/- 0.574938
    a_1 = 3.006151 +/- 0.010033
Sum of squared residuals: 16.755934
Estimated variance of the residuals: 2.094492
Coefficient of determination R^2: 0.999911
F-statistic: 89767.825042
AIC: 15.161673
BIC: 12.069429


The more advanced features
==========================

A lot of interesting information can be extracted from the LinearModel and LinearFit classes.
A simple way to list the available methods on the prompt is:

>>> [method for method in dir(myModel) if not method.startswith("_")]

With the LinearModel class we can compute, for example:

>>> myModel.conditionNumber()
71.994675255533011
>>> myModel.singularValues()
array([ 181.21519026,    2.51706379])
>>> myModel.degreesOfFreedom()
8
>>> myModel.traceHatMatrix()
1.9999999999999996

and more. Some examples with the LinearFit class:

>>> myFit.regressionCoefficients()
array([ 1.51964437,  3.00615141])
>>> myFit.errorBars()
array([ 0.57493781,  0.01003345])
>>> myFit.covarianceMatrix()            # correlation matrix also available
array([[  3.30553480e-01,  -3.49164676e-03],
       [ -3.49164676e-03,   1.00670216e-04]])
>>> myFit.confidenceIntervals(0.05)
(array([ 0.19383541,  2.98301422]), array([ 2.84545333,  3.0292886 ]))

The 95% confidence interval of coefficient a_0 is thus [0.19, 2.85]. See the list 
of methods for more functionality. Note that also weights can be taken into account.


Assessing the models
====================

Test if we can reject the null hypothesis that all coefficients are zero 
with a 5% significance level:

In [57]: myFit.FstatisticTest(0.05)
Out[57]: True

The "True" result means that H0 can indeed be rejected. Now we test if we can
reject the null hypothesis that one particular coefficient is zero with a 5%
significance level

>>> myFit.regressorTtest(0.05)
array([ True,  True], dtype=bool)

We obtained "True" for both coefficients, so both a_0 and a_1 are significant.


Some handy derived classes
==========================

There is a shortcut to define polynomials:

>>> myPolynomial = PolynomialModel(x, "x", [3,1])
>>> print myPolynomial
Model: y = a_0 * x^3 + a_1 * x^1
Expected number of observations: 10

and to define harmonic models:

>>> myHarmonicModel = HarmonicModel(x, "x", [2.13, 3.88], ["f1", "f2"], maxNharmonics=2)
>>> print myHarmonicModel
Model: y = a_0 * sin(2*pi*1*f1*x) + a_1 * cos(2*pi*1*f1*x) + a_2 * sin(2*pi*2*f1*x) + a_3 * cos(2*pi*2*f1*x) + a_4 * sin(2*pi*1*f2*x) + a_5 * cos(2*pi*1*f2*x) + a_6 * sin(2*pi*2*f2*x) + a_7 * cos(2*pi*2*f2*x)
Expected number of observations: 10

One can even add the two models:

>>> myCombinedModel = myPolynomial + myHarmonicModel
>>> print myCombinedModel
Model: y = a_0 * x^3 + a_1 * x^1 + a_2 * sin(2*pi*1*f1*x) + a_3 * cos(2*pi*1*f1*x) + a_4 * sin(2*pi*2*f1*x) + a_5 * cos(2*pi*2*f1*x) + a_6 * sin(2*pi*1*f2*x) + a_7 * cos(2*pi*1*f2*x) + a_8 * sin(2*pi*2*f2*x) + a_9 * cos(2*pi*2*f2*x)
Expected number of observations: 10


Notes
=====
The package is aimed to be robust and efficient.
    - B{Robust}: Singular value decomposition is used to do the fit, where the
                 singular values are treated to avoid numerical problems. This is
                 especially useful when you tend to overfit the data.
    - B{Efficient}: asking twice for the same information will not lead to a
                    recomputation. Results are kept internally. The reason for
                    having a separate LinearModel and LinearFit class is that
                    the bulk of the fitting computations can actually be done
                    without the observations. Fitting the same model to different
                    sets of observations can therefore be done a lot faster this
                    way.
              
"""


import sys
import copy
from math import sqrt,log,pi
import numpy as np
import scipy as sp
import scipy.stats as stats





class LinearModel(object):

    """
    A linear model class 
    
    The class implements a linear model of the form
    y(x) = a_0 * f_0(x) + a_1 * f_1(x) + ... + a_{n-1} f_{n-1}(x)
        - y are the responses (observables)
        - x are the covariates
        - a_i are the regression coefficients (to be determined)
        - f_i(x) are the regressors evaluated in the covariates.
             
    Note: LinearModel instances can be added and printed.
    """


    def __init__(self, regressors, nameList, covMatrix = None, regressorsAreWeighted = False):
    
        """
        Initialiser of the LinearModel class. 
        
        @param regressors: either a list of equally-sized numpy arrays with the regressors 
                           evaluated in the covariates: [f_0(x),f_1(x),f_2(x),...],
                           or an N x M design matrix (numpy array) where these regressor arrays 
                           are column-stacked, with N the number of regressors, and M the number
                           of data points.
        @type regressors: list or ndarray
        @param nameList: list of strings of the regressors names. E.g. ["sin(x)", "cos(x)"]
        @type nameList: list
        @param covMatrix: square array containing the covariance matrix of the
                          observations. Weights, and correlations can hence be
                          taken into account.
        @type covMatrix: ndarray
        @param regressorsAreWeighted: False if regressors are not yet weighted, True otherwise
        @type regressorsAreWeighted: boolean
        @return: a LinearModel instance
        @rtype: LinearModel
        
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([sin(x),x*x], ["sin(x)", "x^2"])
        
        """

        # Sanity check of the 'regressors'
        
        if isinstance(regressors, list):
            self._designMatrix = np.column_stack(np.double(regressors))
            self._nParameters = len(regressors)
            self._nObservations = len(regressors[0])
        elif isinstance(regressors, np.ndarray):
            self._designMatrix = np.double(regressors)
            self._nParameters = regressors.shape[1]
            self._nObservations = regressors.shape[0]
        else:
            raise TypeError, "LinearModel only accepts a list of regressors, or a design matrix"


        # Sanity check of the 'nameList'
    
        if len(nameList) != self._nParameters:
            raise ValueError, "Number of names not equal to number of regressors"
        else:
            self._regressorNames = copy.copy(nameList)

 
        # Sanity check of the 'covMatrix'
        
        if covMatrix == None:
            self._covMatrixObserv = None
        elif not isinstance(covMatrix, np.ndarray):
            raise TypeError, "Covariance matrix of observations needs to be an ndarray object"
        else:
            if len(covMatrix.shape) != 2:
                raise TypeError, "Covariance matrix not a 2-dimensional array"
            elif covMatrix.shape[0] != self._nObservations:
                raise TypeError, "Size of covariance matrix not compatible with number of observations"
            elif covMatrix.shape[1] != covMatrix.shape[0]:
                raise TypeError, "Covariance matris is not a square matrix"
            else:
                self._covMatrixObserv = covMatrix

            
        # If a covariance matrix is specified, and if not already done, compute the weighted 
        # design matrix. The weights are determined by the cholesky decomposition of the 
        # inverse of the covariance matrix. A simple way would be to first invert the covariance 
        # matrix, cholesky-decompose the result, and then do a matrix multiplication. A faster 
        # method is to cholesky-decompose the covariance matrix, and then determine the
        # weighted design matrix by a standard solve.
        # 'choleskyLower' is the lower triangular matrix of the cholesky decomposition.
        
        if (self._covMatrixObserv != None) and (regressorsAreWeighted == False):
            self._choleskyLower = np.linalg.cholesky(self._covMatrixObserv)
            for n in range(self._designMatrix.shape[1]):
                self._designMatrix[:,n] = np.linalg.solve(self._choleskyLower, self._designMatrix[:,n])
        else:
            self._choleskyLower = None
            
        
        
        # Initialisation of several internal variables
        
        self._pseudoInverse = None
        self._standardCovarianceMatrix = None
        self._conditionNumber = None
        self._hatMatrix = None
        self._traceHatMatrix = None
        self._degreesOfFreedom = None
        self._withIntercept = None
        self._U = None                        # svd: A = U * w * V^t
        self._V = None
        self._w = None
        self._B = None
        self._svdTOL = np.finfo(np.double).eps * self._nObservations       # tolerance to ignore singular values
       
        



    def __add__(self, linearModel):
    
        """
        Defines operator '+' of two LinearModel instances.
                
        @param linearModel: a LinearModel instance with the same number of observations as the
                            linear model in 'self'
        @type linearModel: LinearModel class
        @return: a newly created LinearModel object representing a linear model that is the sum 
                 of 'self' and 'linearModel'. 'self' and 'linearModel' are unaltered.
        @rtype: LinearModel
        
        Example:
        
        >>> x = np.linspace(0,10,100)
        >>> lm1 = LinearModel([x], ["x"])
        >>> lm2 = LinearModel([x**2], ["x^2"])
        >>> lm3 = lm1 + lm2    # the same as LinearModel([x, x**2], ["x", "x^2"])
        
        """
        
        if not isinstance(linearModel, LinearModel):
            raise TypeError, "Only a LinearModel can be added to a LinearModel"

        if linearModel.designMatrix().shape[0] != self._nObservations:
            raise ValueError, "Linear model has incompatible design matrix"

        if not np.alltrue(linearModel._covMatrixObserv == self._covMatrixObserv):
            raise ValueError, "Linear model has a different covariance matrix"
            
        designMatrix = np.hstack([self._designMatrix, linearModel.designMatrix()])
        regressorNames = self._regressorNames + linearModel.regressorNames()
        
        if self._covMatrixObserv == None:
            return LinearModel(designMatrix, regressorNames)
        else:
            return LinearModel(designMatrix, regressorNames, self._covMatrixObserv, regressorsAreWeighted = True)
        
        

 
        
               
                      
                             

        

    def designMatrix(self):
    
        """
        Returns the design matrix of the linear model
        
        @return: a N x M design matrix. If the model is written as  
                 y(x) = a_0 * f_0(x) + a_1 * f_1(x) + ... + a_{n-1} f_{n-1}(x)
                 then the design matrix A is defined by A_{i,j} = f_i(x_j). N is
                 the number of regressors, M is the number of data points.
        @rtype: ndarray
                 
        Example:
        
        >>> x = np.linspace(0,10,5)
        >>> lm = LinearModel([x], ["x"])
        >>> lm.designMatrix()
        array([[   0.  ,    0.  ],
        [   2.5 ,    6.25],
        [   5.  ,   25.  ],
        [   7.5 ,   56.25],
        [  10.  ,  100.  ]])

        """
        
        return self._designMatrix












    def _singularValueDecompose(self):
    
        """
        Performs and stores the singular value decomposition of the design matrix
        
        Private class method, not to be used by the user.
        
        @return: nothing
        
        The singular value decomposition of the design matrix is defined as 
             X = U * W * V^t 
        where 
            - X: M x N design matrix.
            - U: M x M unitary matrix 
            - W: M x N diagonal matrix
            - V: N x N unitary matrix
        with M the number of observations, and N the number of regressors. Note that for a 
        unitary matrix A: A^{-1} = A^t
        """
        
        # Compute the singular value decomposition of the design matrix.
        # As the number of singular values is between 1 and min(M,N), only compute
        # the first min(M,N) columns of U, as only these columns are ever used
        # (i.e. multiplied with non-zero values). In the following decomposition:
        # - w is a vector containing the non-zero diagonal elements of W
        # - V^t is returned rather than V
             
        
        self._U,self._w, Vt = np.linalg.svd(self._designMatrix, full_matrices=0)
        self._V = Vt.T    # svd() actually returns V^t rather than V

        # Compute the inverse of the singular values. But set the inverse value
        # to zero if the 
        
        wmax = self._w.max()
        self._invSingularValues = np.zeros_like(self._w)
        for n in range(len(self._w)):
            if self._w[n] > self._svdTOL * wmax:
                self._invSingularValues[n] = 1.0 / self._w[n]


   







    def pseudoInverse(self):
    
        """
        Computes the pseudo-inverse
        
        @return: A numerical safe version of the N x M matrix 
                 X^t X)^{-1} X^t where X is the M x N design matrix, 
                 ^t the transpose, ^{-1} the inverse matrix. 
                 The multiplications are matrix multiplications.
        @rtype: ndarray
                 
        Remarks:
            - The "safeness" is obtained by using a singular value decomposition 
              of X, and treating the singular values in order to avoid numerical 
              problems.
        
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.pseudoInverse()
        array([[  8.95460520e-17,   1.63870968e-01,   1.98709677e-01,
                  1.04516129e-01,  -1.18709677e-01],
               [ -1.07455262e-17,  -1.80645161e-02,  -2.06451613e-02,
                 -7.74193548e-03,   2.06451613e-02]])

        """
        
        if self._pseudoInverse is not None:
            return self._pseudoInverse
        else:
            if self._U == None: self._singularValueDecompose()
            M = self._U.shape[0]      # number of observations
            N = self._U.shape[1]      # number of regressors
            temp = np.zeros((N,M))
            for n in range(N):
                temp[n] = self._invSingularValues[n] * self._U[:,n]
            self._pseudoInverse = np.dot(self._V, temp)
            
            return self._pseudoInverse













    def standardCovarianceMatrix(self):
    
        """
        Computes the standard covariance matrix.
                 
        @return: A  numerical safe version of C = (X^t X)^{-1} where
                 X is the design matrix, ^t the transpose, ^{-1} the 
                 inverse matrix. The multiplications are matrix 
                 multiplications. 
        @rtype: ndarray
         
        Remarks:
            - The "safeness" is obtained by using a singular value 
              decomposition of X, and treating the singular values to avoid 
              numerical problems.
            - The matrix C is the covariance matrix of the regression 
              coefficients if the signal noise would be i.i.d. with sigma = 1.
            - In the case of sigma != 1, the actual covariance matrix
              can be obtained by simply rescaling the standard covariance
              matrix. See L{covarianceMatrix}.
                 
        Example: 
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.standardCovarianceMatrix()
        array([[ 0.09135484, -0.01032258],
               [-0.01032258,  0.00123871]])

        """
        
        if self._standardCovarianceMatrix != None:
            return self._standardCovarianceMatrix
        else:
            if self._V == None: self._singularValueDecompose()
            self._standardCovarianceMatrix = np.dot(self._V, np.dot(np.diag(self._invSingularValues**2), self._V.T))
            return self._standardCovarianceMatrix
        










    def conditionNumber(self):
    
        """
        Computes the condition number.
                         
        @return: The ratio of the highest to the lowest singular value,
                 where the singular values are obtained by singular value 
                 decomposition of the design matrix. A large condition number 
                 implies an ill-conditioned design matrix.
        @rtype: double
        
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm..conditionNumber()
        35.996606504814814
        
        """
        
        if self._conditionNumber is not None:
            return self._conditionNumber
        else:
            if self._w == None: self._singularValueDecompose()
            self._conditionNumber = self._w.max() / self._w.min()
            return self._conditionNumber
        










    def singularValues(self):
    
        """
        Return the non-zero singular values 
        as obtained with the singular value decomposition
        
        @return: The non-zero singular values.
        @rtype: ndarray
       
        Remarks: 
            - the singular values can be used a diagnostic to see how 
              ill-conditioned the matrix inversion is.
              
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.singularValues()
        array([ 118.34194851,    3.28758625])

        """
        
        if self._w != None:
            return self._w
        else:
            self._singularValueDecompose()
            return self._w












    def hatMatrix(self):
    
        """
        Computes the hat matrix.
        
        The hat matrix is defined by H = X (X^t X)^{-1} X^t with
            - X the (weighted) design matrix
            - ^t the transpose, 
            - ^{-1} the inverse matrix 
        and where all multiplications are matrix multiplications. 
        It has the property that \hat{y} = H * y where
            - \hat{y} are the (weighted) predictions
            - y the (weighted) observations. 
        So, it "puts a hat" on the (weighted) observation vector.
        
        @return: The hat matrix.
        @rtype: ndarray
        
        Remarks:
            - as the number of data points may be large, this matrix may
              become so huge that it no longer fits into your computer's
              memory
            - for weighted (and/or correlated) observations, the resulting 
              weighted hat matrix does not put a hat on the unweighted 
              observations vector, only on the weighted one.
              
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.hatMatrix()
        array([[  9.32150093e-32,   1.56705591e-16,   1.79092104e-16,
                  6.71595390e-17,  -1.79092104e-16],
               [  1.56705591e-16,   2.96774194e-01,   3.67741935e-01,
                  2.12903226e-01,  -1.67741935e-01],
               [  1.79092104e-16,   3.67741935e-01,   4.77419355e-01,
                  3.29032258e-01,  -7.74193548e-02],
               [  6.71595390e-17,   2.12903226e-01,   3.29032258e-01,
                  3.48387097e-01,   2.70967742e-01],
               [ -1.79092104e-16,  -1.67741935e-01,  -7.74193548e-02,
                  2.70967742e-01,   8.77419355e-01]])
                  
        """
        
        if self._hatMatrix is not None:
            return self._hatMatrix
        else:
            if self._U is None: self._singularValueDecompose()
            self._hatMatrix = np.dot(self._U, self._U.T)
            return self._hatMatrix
        











    def traceHatMatrix(self):
    
        """
        Returns the trace of the (possibly weighted) hat matrix.
        
        The computation of the entire MxM hat matrix (M the number of 
        observations) is avoided to calculate the trace. This is useful 
        when the number of observations is so high that the hat matrix 
        does no longer fit into the memory.
        
        @return: The trace of the hat matrix.
        @rtype: ndarray
        
        Remarks:
            - If the hat matrix was computed anyway it will be used.
            - The trace of the hat matrix that puts a hat on the original
              observations equals the trace of the hat matrix that puts a 
              hat on the weighted observations.
              
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.hatMatrix()
        1.9999999999999993
        
        """
        
        if self._traceHatMatrix is not None:
            return self._traceHatMatrix
        else:
            if self._hatMatrix is not None:
                self._traceHatMatrix = self._hatMatrix.trace()
            else:
                if self._U is None: self._singularValueDecompose()
                self._traceHatMatrix = 0.0
                for k in range(self._U.shape[1]):
                    self._traceHatMatrix += np.sum(np.square(self._U[:,k]))
            return self._traceHatMatrix
            
        











    def degreesOfFreedom(self):
    
        """
        Returns the effective number of degrees of freedom.
        
        The d.o.f. is defined as the trace of the hat matrix. See also:
        U{Wikipedia<http://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)>}
                 
        @return: The effective number of degrees of freedom.
        @rtype: integer
        
        Example:
        
        >>> x = linspace(0,10,5)
        >>> lm = LinearModel([x, x**2], ["x", "x^2"])
        >>> lm.degreesOfFreedom()
        3
        
        """
        
        if self._degreesOfFreedom is not None:
           return self._degreesOfFreedom
        else:
           self._degreesOfFreedom = self._nObservations - int(round(self.traceHatMatrix()))
           return self._degreesOfFreedom
                   










    def regressorNames(self):
    
        """
        Returns the regressor names.
        
        @return: list with strings containing the names of the regressors.
        @rtype: list
    
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([sin(x),x*x], ["sin(x)", "x^2"]) 
        >>> lm.regressorNames()
        ["sin(x)", "x^2"]
        
        """
        
        return self._regressorNames
    
    












    def containsRegressorName(self, regressorName):
    
        """
        Checks if the linear model contains the given regressor.
        
        @param regressorName: name of a regressor
        @type regressorName: string
        @return: True if one of the regressor names is identical 
                 to the input 'regressorName', else False.
        @rtype: boolean
        
        Example: 
           
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([sin(x), x*x], ["sin(x)", "x^2"])
        >>> lm.containsRegressorName("x^2")
        True
        >>> lm.containsRegressorName("cos(x)")
        False
        
        """
        
        return regressorName in self._regressorNames
        
        







    def someRegressorNameContains(self, someString):
    
        """
        Checks if a regressor name contains a given substring.
        
        @param someString: substring to be checked
        @type someString: string
        
        @return: True if at least one of the regressor names contains
                 the substring someString, else return false.
        @rtype: boolean
                         
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([sin(x), x*x], ["sin(x)", "x^2"])
        >>> lm.someRegressorNameContains("sin")
        True
        
        """

        for name in self._regressorNames:
            if someString in name: return True
            
        return False










    def withIntercept(self):
    
        """
        Checks if there is a constant regressor
        
        @return: True if there is a constant regressor (the intercept).
                 False otherwise.
        @rtype: boolean        
                
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([1, x*x], ["1", "x^2"])
        >>> lm.withIntercept()
        True
        
        """
        
        if self._withIntercept is not None:
            return self._withIntercept
        else:
            self._withIntercept = False
            for n in range(self._nParameters):
                if len(np.unique(self._designMatrix[:,n])) == 1:
                    self._withIntercept = True
                    break
            return self._withIntercept

    








    def nObservations(self):
    
        """
        Returns the number of observations
        
        @return: the number of observations that the LinearModel expects
        @rtype: integer
        
        Example:
        
        >>> x = linspace(0,10,23)
        >>> lm = LinearModel([x], ["x"])
        >>> lm.nObservations()
        23
        
        """
        
        return self._nObservations
        
        
 
        
        







    def nParameters(self):
    
        """
        Returns the number of fit coefficients
        
        @return: the numbero fit coefficients of the linear model
        @rtype: integer
        
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([sin(x),x*x], ["sin(x)", "x^2"])
        >>> lm.nParameters()
        2
          
        """
        
        return self._nParameters
        
        










    def fitData(self, observations):
        
        """
        Create and return a LinearFit object.
        
        From this object the user can extract all kinds of quantities 
        related to the linear fit. See L{LinearFit}.
        
        @param observations: array with the observations y_i.  The size 
                             of the array should be compatible with the number
                             of expected observations in the LinearModel.
        @type observations: ndarray
        @return: a LinearFit instance, containing the fit of the observations
                 with the linear model.
        @rtype: LinearFit
                                                        
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([x], ["x"])
        >>> obs = x + normal(0.0, 1.0, 100)         # some simulated observations
        >>> linearfit = lm.fitData(obs)
        
        """
        
        if len(observations) != self._nObservations:
            raise ValueError, "Number of observations should be %d != %d" % (self._nObservations, len(observations))
            
        return LinearFit(self, observations)
        
        









    
    
    def submodels(self, Nmin = 1, Nmax = None, nested = True, ranks = None, simpleFirst=True):
    
        """
        Returns a generator capable of generating submodels
        
        The generator sequentially returns the submodels which are instances
        of LinearModel with only a subset of the regressors of the parent model.
        The submodels with fewer regressors will always appear first in the list. 
        The order in which the regressors appear in the submodels will be from 
        low rank to high rank.
        
        @param Nmin: Minimum number of regressors (or groups of equally ranking 
                     regressors). Should be minimal 1. 
        @type Nmin: integer
        @param Nmax: Maximum number of regressors (or groups of equally ranking 
                     regressors). If == None, it will take the maximum number 
                     possible. Note: if the ranks = [0,1,2,2] then there are 
                     maximum 3 "groups" of regressors, as the latter 2 regressors 
                     will always be taken together.
        @type Nmax: integer
        @param nested: if False: return all possible submodels, i.e. all possible 
                       combinations of regressors. If true: return only nested
                       submodels with combinations of regressors (or groups of 
                       equally ranked regressors) so that a higher ranked regressors 
                       will never appear without the lower ranked regressors
        @type nested: boolean
        @param ranks: list with (integer) rankings of the regressors, only the 
                      relative ranking of the numbers will be used, so the exact 
                      numbers are not important. For example, [1,2,3] has the 
                      same result as [10,24,35]. Regressors having the same rank 
                      should always be included together in a submodel, regardless 
                      how the input parameter 'nested' has been set. In case of 
                      nested submodels, a regressor with a higher rank should 
                      never occur without the regressors with lower rank in the 
                      same submodel.
        @type ranks: list with integers
        @param simpleFirst: if True: generate the more simple submodels first.
                            if False: generate the complicated submodels first.
        @type simpleFirst: boolean
        @return: a python generator implementing the yield method
        @rtype: generator
        
        Examples: 
        If the regressors are [1, x, x**2] then
            - Nmin=1, nested=False, ranks=None gives the submodels
              [1], [x], [x**2], [1,x], [1,x**2], [x,x**2], [1,x,x**2] 
            - Nmin=1, Nmax=2, nested=False, ranks=None gives the submodels
              [1], [x], [x**2], [1,x], [1,x**2], [x,x**2]
            - Nmin=2, nested=False, ranks=None gives the submodels
              [1,x], [1,x**2], [x,x**2], [1,x,x**2]
            - Nmin=3, nested=True, ranks=None gives the submodel
              [1,x,x**2]
            - Nmin=3, nested=True, ranks=[1,0,3] gives the submodel
              [x,1,x**2]  
            - Nmin=2, nested=False, ranks=[1,0,3] gives the submodels
              [x,1], [x,x**2], [1,x**2], [x,1,x**2]
            - Nmin=2, nested=True, ranks=[1,0,3] gives the submodels
              [x,1], [x,1,x**2]
        If the regressors are [1,x,sin(2x),cos(2x)] then
            - Nmin=3, nested=True, ranks=[0,1,2,2] gives the submodel
              [1,x,sin(2x),cos(2x)]
            - Nmin=2, nested=True, ranks=[0,1,2,2] gives the submodels
              [1,x], [1,x,sin(2x),cos(2x)]
            - Nmin=2, nested=False, ranks=[0,1,2,2] gives the submodels
              [1,x], [1,sin(2x),cos(2x)], [x,sin(2x),cos(2x)], [1,x,sin(2x),cos(2x)]
            - Nmin=1, Nmax=2, nested=False, ranks=[0,1,2,2] gives the submodels
              [1], [x], [sin(2x),cos(2x)], [1,x], [1,sin(2x),cos(2x)], [x,sin(2x),cos(2x)]
              
        >>> 
        """
        
        # If no ranking of the regressors are given, consider the first regressor
        # having the most important, and going down, the last regressor the least 
        # important
        
        if ranks == None:
            ranks = np.arange(self._nParameters)
            
        # Sort the (unique) ranks from low (most important) to high (least important)    
        
        uniqueRanks = np.unique(ranks)
        if Nmax == None:
            Nmax = len(uniqueRanks)
        
        # nRegressor is a list to be looped over, containing for each submodel the number
        # of regressors. Be aware that a group of regressors with equal rank are considered to
        # be only 1 regressor in this accounting. 
        
        if simpleFirst:
            nRegressorRange = range(Nmin,Nmax+1)
        else:
            nRegressorRange = range(Nmax,Nmin-1,-1)
        
        # Make a generator yield the submodels
        
        if nested == True:
            for n in nRegressorRange:
                nameListSub = []
                regressorListSub = []
                indices = [k for m in uniqueRanks[:n] for k in np.where(ranks == m)[0]] 
                nameListSub += [self._regressorNames[k] for k in indices]
                regressorListSub += [self._designMatrix[:,k] for k in indices]
                if self._covMatrixObserv != None:
                    covMatrixObservSub = self._covMatrixObserv[indices,:][:,indices]
                    yield LinearModel(regressorListSub, nameListSub, covMatrixObservSub, regressorsAreWeighted=True)
                else:
                    yield LinearModel(regressorListSub, nameListSub)
        else:
            comboList = [combo for comboSize in nRegressorRange for combo in combinations(uniqueRanks, comboSize)]
            for combo in comboList:
                indices = [k for n in combo for k in np.where(ranks==n)[0]]
                nameListSub = [self._regressorNames[k] for k in indices]
                regressorListSub = [self._designMatrix[:,k] for k in indices]
                if self._covMatrixObserv != None:
                    covMatrixObservSub = self._covMatrixObserv[indices,:][:,indices]
                    yield LinearModel(regressorListSub, nameListSub, covMatrixObservSub, regressorsAreWeighted=True)
                else:
                    yield LinearModel(regressorListSub, nameListSub)
                
    







    

    def copy(self):
        
        """
        Returns a duplicate of the current LinearModel
        
        @return: of duplicate of self. Not a reference,
                 a real copy.
        @rtype: LinearModel
                 
        Example:
        
        >>> x = linspace(0,10,100)
        >>> lm = LinearModel([x], ["x"])
        >>> lm2 = lm.copy()
        
        """
        
        if self._covMatrixObserv != None:
            return LinearModel(self._designMatrix.copy(), copy.copy(self._regressorNames), self._covMatrixObserv.copy())
        else:
            return LinearModel(self._designMatrix.copy(), copy.copy(self._regressorNames))
            

        
        





    
    
     
    def __str__(self):
    
        """
        Generates what is written to the console if the LinearModel is printed

        @return: nothing
        
        """
        
        if self._regressorNames[0] == "1":
            str = "Model: y = a_0"
        else:
            str = "Model: y = a_0 * %s" % self._regressorNames[0]
        
        for n in range(1, self._nParameters):
            if self._regressorNames[n] == "1":
                str += " + a_%d" % n
            else:
                str += " + a_%d * %s" % (n, self._regressorNames[n])
                
        str += "\nExpected number of observations: %d" % self._nObservations
        
        return str
        
        












class PolynomialModel(LinearModel):

    """
    Class to implement a polynomial model
    
    The class implements a polynomial model function of the form
        y(x) = a_0 * a_1 * x + a_2 * x^2 + ...
    where y are the responses (observables), x are the covariates, 
    the a_i are the regression coefficients (to be determined).
    Such a class provides convenient way of creating a polynomial model.
             
    Remarks:
        - PolynomialModel instances can be added to other PolynomialModel
          or LinearModel instances to give a new LinearModel instance.
    """
    
    
    def __init__(self, covariate, covariateName, orderList, covMatrix = None):

        """
        Initialiser of the PolynomialModel class. 
        
        @param covariate: an array of size N with the covariate values, 
                          where N is the number of observations
        @type covariate: ndarray
        @param covariateName: the (preferably short) name of the covariate. 
                              E.g. "t", or "x".
        @type covariateName: string
        @param orderList: list of the orders that should be included in the
                          model. E.g. [3,0,1] implies the model
                          a_0 * x^3 + a_1 + a_2 * x 
        @type orderList: list
        @param covMatrix: square array containing the covariance matrix of the
                          observations. This way, weights and correlations can 
                          be taken into account.
        @type covMatrix: ndarray
        @return: an instance of PolynomialModel, subclass of LinearModel
        @rtype: PolynomialModel
               
        Example:
        
        >>> time = linspace(0,10,100)
        >>> pm = PolynomialModel(time, "t", [2,1])
        
        """
        
        # Create and fill the design matrix
        
        designMatrix = np.empty((len(covariate), len(orderList)))
        for n in range(len(orderList)):
            designMatrix[:,n] = covariate**orderList[n]
 
               
        # List the regressor names
        
        regressorNames = []
        for n in range(len(orderList)):
            if orderList[n] == 0:
                regressorNames += ["1"]
            else:
                regressorNames += ["%s^%d" % (covariateName, orderList[n])]
 
 
        # Let the base class do the initialising
 
        super(PolynomialModel, self).__init__(designMatrix, regressorNames, covMatrix, regressorsAreWeighted=False)
                
                
        # Testing for an intercept can be done more efficiently for polynomials
        
        if 0 in orderList:
            self._withIntercept = True
        else:
            self._withIntercept = False
            




class HarmonicModel(LinearModel):

    """
    Class to implement a harmonic model function.
    
    The class implements a harmonic model of the form
        y(x) =   a_0 * sin(2*pi*f1*t) + a_1 * cos(2*pi*f1*t)
               + a_2 * sin(2*pi*f2*t) + a_3 * cos(2*pi*f2*t) + ...
    where y are the responses (observables), x are the covariates, 
    the a_i are the regression coefficients (to be determined).
    Such a class provides a more convenient short way of creating a
    harmonic model.
    
    Remarks: 
        - HarmonicModel instances can be added to other HarmonicModel
          or LinearModel instances to give a new LinearModel instance.
    """
    
    
    def __init__(self, covariate, covariateName, freqList, freqNameList, 
                       maxNharmonics=1, covMatrix = None):

        """
        Initialiser of HarmonicModel class, with *known* frequencies. 
        
        @param covariate: an array of size N with the covariate values, 
                          where N is the number of observations.
        @type covariate: ndarray
        @param covariateName: a (preferably short) name of the covariate. 
                              E.g. "t", or "x".
        @type covariateName: string
        @param freqList: the frequencies (unit: inverse unit of covariate)
        @type freqList: list
        @param freqNameList: the (preferably short) names of the frequencies. 
                             E.g. "f_1" or "nu_1".
        @type freqNameList: list
        @param maxNharmonics: integer values > 1. The number of harmonics for
                              each frequency. For example:
                              1 : f_1, f_2, ...
                              2 : f_1, 2*f_1, f_2, 2*f_2, ...
        @type maxNharmonics: list
        @param covMatrix: square array containing the covariance matrix of the
                          observations. Weights, and correlations can hence be
                          taken into account.
        @type covMatrix: ndarray
        @return: an instance of HarmonicModel, subclass of LinearModel
        @rtype: HarmonicModel
             
        Example: 
           
        >>> time = linspace(0,10,100)
        >>> hm = HarmonicModel(time, "t", [2.34, 8.9], ["f_1", "f_2"], 1)
        
        """
        
        # Create and fill the design matrix
        
        nParam = 2*len(freqList)*maxNharmonics
        designMatrix = np.empty((len(covariate), nParam))
        runningIndex = 0
        for n in range(len(freqList)):
            for k in range(1,maxNharmonics+1):
                designMatrix[:,runningIndex] = np.sin(2.*pi*k*freqList[n]*covariate)
                designMatrix[:,runningIndex+1] = np.cos(2.*pi*k*freqList[n]*covariate)
                runningIndex += 2
        
        # List the regressor names
        
        regressorNames = []
        for n in range(len(freqList)):
            for k in range(1,maxNharmonics+1):
                 regressorNames += ["sin(2*pi*%d*%s*%s)" % (k, freqNameList[n], covariateName), 
                                    "cos(2*pi*%d*%s*%s)" % (k, freqNameList[n], covariateName)]
                                             

        # Let the base class do the initialising
 
        super(HarmonicModel, self).__init__(designMatrix, regressorNames, covMatrix, regressorsAreWeighted=False)
                
                
        # There is no intercept for a harmonic model
        
        self._withIntercept = False
















class LinearFit(object):

    """
    A class that allows to easily extract the fit coefficients, 
    covariance matrix, predictions, residuals, etc.
    
    Remarks: LinearFit instances can be printed
    """
    
    
    def __init__(self, linearModel, observations):
    
        """
        Initialises the LinearFit instance
        
        @param linearModel: a LinearModel instance
        @type linearModel: LinearModel
        @param observations: the observations. The size of the array must be 
                             compatible with the number of observations that 
                             the linearModel expects.
        @type observations: ndarray
        @return: a LinearFit instance
        @rtype: LinearFit
        
        """
        
        # Store some ndarrays internally for later
        
        self._linearModel = linearModel 
        self._nObservations = len(observations)
        self._nParameters = linearModel.designMatrix().shape[1]

         
        # Init some quantities as None, signaling that they have not yet been computed
        
        self._regressionCoefficients = None
        self._residuals = None
        self._weightedResiduals = None
        self._predictions = None
        self._weightedPredictions = None
        self._sumSqResiduals = None
        self._sumSqWeightedResiduals = None
        self._covarianceMatrix = None
        self._correlationMatrix = None
        self._errorBars = None
        self._t_values = None
        self._coefficientOfDetermination = None
        self._weightedCoeffOfDetermination = None
        self._residualVariance = None
        self._weightedResidualVariance = None
        self._BICvalue = None
        self._AICvalue = None
        self._Ttest = None
        self._Fstatistic = None
        self._weightedFstatistic = None
        
        
        # If a covariance of the observations was specified, weight the observations.
        
        self._originalObservations = observations
        if linearModel._covMatrixObserv != None:
            self._weightedObservations = np.linalg.solve(linearModel._choleskyLower, observations)
        else:
            self._weightedObservations = observations
        


    def observations(self, weighted=False):
        
        """
        Returns the observations
        
        @param weighted: If false return the original observations, if
                         True, return de decorrelated observations.
        @type weighted: boolean
        @return: the observations that were used to initialize the LinearFit instance
        @rtype: ndarray
                         
        Remarks: 
            - if no covariance matrix of the observations was specified for
              the model, the "decorrelated" observations are identical to the
              original observations.
        """
        
        if weighted:
            return self._weightedObservations
        else:
            return self._originalObservations
        

 
        
    def regressionCoefficients(self):
    
        """
        Returns the regression coefficients. 
        
        @return: the fit coefficients
        @rtype: ndarray
        
        """
    
        if self._regressionCoefficients is not None:
            return self._regressionCoefficients
        else:
            self._regressionCoefficients = np.dot(self._linearModel.pseudoInverse(), 
                                                  self._weightedObservations)
            return self._regressionCoefficients
    
    


    def sumSqResiduals(self, weighted=False):
        
        """
        Returns the sum of the squared residuals
        
        @param weighted: If false the unweighted residuals are used,
                         if True, the weighted ones are used.
        @type weighted: boolean
        @return: the sum of the squared residuals
        @rtype: double
        
        """
        
        if weighted:
            if self._sumSqWeightedResiduals is not None:
                return self._sumSqWeightedResiduals
            else:
                self._sumSqWeightedResiduals = np.sum(np.square(self.residuals(weighted=True)))
                return self._sumSqWeightedResiduals
        else:
            if self._sumSqResiduals is not None:
                return self._sumSqResiduals
            else:
                self._sumSqResiduals = np.sum(np.square(self.residuals(weighted=False)))
                return self._sumSqResiduals

            
   



    def residualVariance(self, weighted=False):
    
        """
        Estimates the variance of the residuals of the time series.
        
        As normalization the degrees of freedom is used
        
        @param weighted: If false the unweighted residuals are used,
                         if True, the weighted ones are used.
        @type weighted: boolean
        @return: the variance of the residuals
        @rtype: double
        
        """
        
        if weighted:
            if self._weightedResidualVariance is not None:
                return self._weightedResidualVariance
            else:
                self._weightedResidualVariance = self.sumSqResiduals(weighted=True)     \
                                                 / self._linearModel.degreesOfFreedom()
                return self._weightedResidualVariance
        else:
            if self._residualVariance is not None:
                return self._residualVariance
            else:
                self._residualVariance = self.sumSqResiduals(weighted=False)            \
                                         / self._linearModel.degreesOfFreedom()
                return self._residualVariance
        



    def covarianceMatrix(self):
    
        """
        Returns the covariance matrix of the fit coefficients
        
        @return: the MxM covariance matrix of the fit coefficients with
                 M the number of fit coefficients
        @rtype: ndarray
        
        """
        
        if self._covarianceMatrix is not None:
            return self._covarianceMatrix
        else:
            self._covarianceMatrix = self._linearModel.standardCovarianceMatrix()       \
                                     * self.residualVariance(weighted=True)
            return self._covarianceMatrix
        
    

    
    def correlationMatrix(self):
    
        """
        Returns the correlation matrix of the fit coefficients
        
        @return: the MxM correlation matrix of the fit coefficients with
                 M the number of fit coefficients
        @rtype: ndarray   
        
        """
        
        if self._correlationMatrix is not None:
            return self._correlationMatrix
        else:
            cov = self.covarianceMatrix()
            self._correlationMatrix = np.empty_like(cov)
            cor = self._correlationMatrix                  # shorter alias
            for n in range(cor.shape[0]):
                for m in range(n):
                    cor[n,m] = cov[n,m] / sqrt(cov[n,n] * cov[m,m])
                    cor[m,n] = cor[n,m]
                cor[n,n] = 1.0
            return self._correlationMatrix
        



    def errorBars(self):
    
        """
        Returns the formal error bars of the fit coefficients
        
        @return: The square root of the diagonal of the covariance matrix
        @rtype: ndarray
        
        """
        
        if self._errorBars is not None:
            return self._errorBars
        else:
            self._errorBars = np.sqrt(self.covarianceMatrix().diagonal())
            return self._errorBars
    


    def confidenceIntervals(self, alpha):
    
        """
        Returns the symmetric (1-alpha) confidence interval around the fit coefficients
        
        E.g. if alpha = 0.05, the 95% symmetric confidence interval is returned.
                                  
        @param alpha: confidence level in ]0,1[
        @type alpha: double
        @return: the arrays lowerBounds[0..K-1] and upperBounds[0..K-1] which contain
                 respectively the lower and the upper bounds of the symmetric interval
                 around the fit coefficients, where K is the number of fit coeffs.
        @rtype: tuple
        
        Remarks:
            - The formula used assumes that the noise on the observations is 
              independently, identically and gaussian distributed.
              
        """

        t = stats.t.ppf(1.0-alpha/2., self._linearModel.degreesOfFreedom())

        lowerBounds = self.regressionCoefficients() - t * self.errorBars()
        upperBounds = self.regressionCoefficients() + t * self.errorBars()
        
        return lowerBounds, upperBounds
        
        
        
    def t_values(self):
    
        """
        Returns the formal t-values of the fit coefficients
        
        @return: t-values = fit coefficients divided by their formal error bars
        @rtype: ndarray
        
        """
        
        if self._t_values is not None:
            return self._t_values
        else:
            self._t_values = self.regressionCoefficients() / self.errorBars()
            return self._t_values
            
    

    
    def regressorTtest(self, alpha):
    
        """
        Performs a hypothesis T-test on each of the regressors
        
        Null hypothesis: H0 : fit coefficient == 0
        Alternative hypothesis : H1 : fit coefficient != 0
        
        @param alpha: significance level of the hypothesis test. In ]0,1[.
                      E.g.: alpha = 0.05
        @type alpha: double
        @return: a boolean array of length K with K the number of regressors. 
                 "True" if the null hypothesis was rejected for the regressor, 
                 "False" otherwise.
        @rtype: ndarray
        
        Remarks: 
            - This test assumes that the noise is independently,
              identically, and gaussian distributed. It is not robust 
              against this assumption.
              
        """
        
        t = stats.t.ppf(1.0-alpha/2., self._linearModel.degreesOfFreedom())
         
        if self._Ttest is not None:
            return self._Ttest
        else:
            self._Ttest = np.fabs(self.t_values()) > t
            return self._Ttest

 

    def predictions(self, weighted=False):
    
        """
        Returns the predicted (fitted) values
        
        It concerns the predictions for the original observations, 
        not the decorrelated ones.
                 
        @param weighted: If True/False, the predictions for the 
                         decorrelated/original observations will be used.
        @type weighted: boolean
                         
        Remarks:
            - If no covariance matrix of the observations was specified 
              for the model, the weighted predictions are identical to 
              the unweighted ones.
        """
        
        if weighted:
            if self._weightedPredictions is not None:
                return self._weightedPredictions
            else:
                if self._linearModel._covMatrixObserv != None:
                    self._weightedPredictions = np.dot(self._linearModel._choleskyLower, np.dot(self._linearModel.designMatrix(), self.regressionCoefficients()))
                else:
                    self._weightedPredictions = np.dot(self._linearModel.designMatrix(), self.regressionCoefficients())
                return self._weightedPredictions
        else:
            if self._predictions is not None:
                return self._predictions
            else:
                if self._linearModel._covMatrixObserv != None:
                    self._predictions = np.dot(self._linearModel._choleskyLower, np.dot(self._linearModel.designMatrix(), self.regressionCoefficients()))
                else:
                    self._predictions = np.dot(self._linearModel.designMatrix(), self.regressionCoefficients())
                return self._predictions

    


    def residuals(self, weighted=False):
    
        """
        Returns an array with the residuals.
        
        Residuals = observations minus the predictions
        
        @param weighted: If True/False, the residuals for the 
                         decorrelated/original observations will be used.
        @type weighted: boolean
                         
        Remarks:
            - If no covariance matrix of the observations was specified for the
              model, the weighted residuals are identical to the unweighted ones.
              
        """
        
        if weighted:
            if self._weightedResiduals is not None:
                return self._weightedResiduals
            else:
                self._weightedResiduals = self._weightedObservations - self.predictions(weighted=True)
                return self._weightedResiduals
        else:
            if self._residuals is not None:
                return self._residuals
            else:
                self._residuals = self._originalObservations - self.predictions(weighted=False)
                return self._residuals

            
     


    def coefficientOfDetermination(self, weighted=False):
    
        """
        Returns the coefficient of determination
        
        The coeff of determination is defined by  1 - S1 / S2 with
            - S1 the sum of squared residuals:
              S1 = \frac{\sum_{i=1}^N (y_i - \hat{y}_i)^2
            - S2 the sample variance w.r.t. the mean:
              S2 = \frac{\sum_{i=1}^N (y_i - \overline{y})^2
            - If there is no intercept term in the model, then S2 is computed
              S2 = \frac{\sum_{i=1}^N (y_i)^2
        
        @param weighted: If True/False, the residuals for the 
                         decorrelated/original observations will be used.
        @type weighted: boolean
        @return: coefficient of determination
        @rtype: double
        
        Remarks:
            - If no covariance matrix of the observations was specified 
              for the model, the weighted coeff of determination is 
              identical to the unweighted one.
        """
        
        if weighted:
            if self._weightedCoeffOfDetermination is not None:
                return self._weightedCoeffOfDetermination
            else:
                if self._linearModel.withIntercept():
                    sampleVariance = np.sum(np.square(self._weightedObservations - np.mean(self._weightedObservations)))
                    self._weightedCoeffOfDetermination = 1.0 - self.sumSqResiduals(weighted=True) / sampleVariance
                else:
                    sampleVariance = np.sum(np.square(self._weightedObservations))
                    self._weightedCoeffOfDetermination = 1.0 - self.sumSqResiduals(weighted=True) / sampleVariance
                
                return self._coefficientOfDetermination
        else:
            if self._coefficientOfDetermination is not None:
                return self._coefficientOfDetermination
            else:
                if self._linearModel.withIntercept():
                    sampleVariance = np.sum(np.square(self._originalObservations - np.mean(self._originalObservations)))
                    self._coefficientOfDetermination = 1.0 - self.sumSqResiduals(weighted=False) / sampleVariance
                else:
                    sampleVariance = np.sum(np.square(self._originalObservations))
                    self._coefficientOfDetermination = 1.0 - self.sumSqResiduals(weighted=False) / sampleVariance
                
                return self._coefficientOfDetermination

                



    def BICvalue(self):
    
        """
        Returns the Bayesian Information Criterion value. 
        
        @return: BIC value
        @rtype: double
        
        Remarks:
            - Also called the Schwartz Information Criterion (SIC)
            - Gaussian noise is assumed, with unknown variance sigma^2
            - Constant terms in the log-likelihood are omitted
            
        TODO: . make a weighted version
        """

        if self._BICvalue is not None:
            return self._BICvalue
        else:
            # Include the sigma of the noise in the number of unknown fit parameters
            nParam = self._nParameters + 1  
            self._BICvalue = self._nObservations * log(self.sumSqResiduals(weighted=False)/self._nObservations)    \
                             + nParam * log(self._nObservations)
            return self._BICvalue
            
    
    
    
    def AICvalue(self):
    
        """
        Returns the 2nd order Akaike Information Criterion value 
        
        @return: AIC value
        @rtype: double
        
        Remarks:
            - Gaussian noise is assumed, with unknown variance sigma^2
            - Constant terms in the log-likelihood were omitted
            
        TODO: . make a weighted version
        """

        if self._AICvalue is not None:
            return self._AICvalue
        else:
            # Include the sigma of the noise in the number of unknown fit parameters
            nParam = self._nParameters + 1
            self._AICvalue = self._nObservations * log(self.sumSqResiduals(weighted=False)/self._nObservations)     \
                       + 2*nParam + 2*nParam*(nParam+1)/(self._nObservations - nParam - 1)
            return self._AICvalue
        



    def Fstatistic(self, weighted=False):
    
        """
        Returns the F-statistic, commonly used to assess the fit.
        
        @param weighted: If True/False, the residuals for the 
                         decorrelated/original observations will be used.
        @type weighted: boolean
        @return: the F-statistic
        @rtype: double
                         
        Remarks:
            - if no covariance matrix of the observations was specified for the
              model, the weighted F-statistic is identical to the unweighted one

        """
        
        if weighted:
            if self._weightedFstatistic is not None:
                return self._weightedFstatistic
            else:
                Rsq = self.coefficientOfDetermination(weighted=True)
                df = self._linearModel.degreesOfFreedom()
            
                if self._nParameters == 1:
                    self._Fstatistic = Rsq / (1-Rsq) * df / self._nParameters
                else:
                    self._Fstatistic = Rsq / (1-Rsq) * df / (self._nParameters-1)
            
                return self._weightedFstatistic            
        else:
             if self._Fstatistic is not None:
                return self._Fstatistic
             else:
                Rsq = self.coefficientOfDetermination(weighted=False)
                df = self._linearModel.degreesOfFreedom()
            
                if self._nParameters == 1:
                    self._Fstatistic = Rsq / (1-Rsq) * df / self._nParameters
                else:
                    self._Fstatistic = Rsq / (1-Rsq) * df / (self._nParameters-1)
            
                return self._Fstatistic                       
    
    

    def FstatisticTest(self, alpha, weighted=False):
    
        """
        Performs a hypothesis F-test on the fit
        
        Null hypothesis: H0 : all fit coefficients == 0
        Alternative hypothesis : H1 : at least one fit coefficient != 0
                 
        Stated otherwise (with R^2 the coefficient of determination):
        Null hypothesis: H0 : R^2 == 0
        Alternative hypothesis: H1: R^2 != 0
                 
        @param alpha: significance level of the hypothesis test. In ]0,1[.
        @type alpha: double
        @param weighted: If True/False, the weighted/not-weighted 
                         F-statistic will be used.
        @type weighted: boolean
        @return: "True" if null hypothesis was rejected, "False" otherwise
        @rtype: boolean
                            
        Remarks:
            - if no covariance matrix of the observations was specified for the
              model, the weighted F-test is identical to the unweighted one
              
        """
        
        df = self._linearModel.degreesOfFreedom()
        
        if self._nParameters == 1:
            Fcritical = stats.f.ppf(1.0-alpha, self._nParameters, df)
        else:
            Fcritical = stats.f.ppf(1.0-alpha, self._nParameters-1, df)

        return (self.Fstatistic(weighted) > Fcritical)



    def summary(self, outputStream = sys.stdout):
    
        """
        Writes some basic results of fitting the model to the observations.
        
        @param outputStream: defaulted to the standard output console, but can be
                             replaced by another open stream, like a file stream.
        @type outputStream: stream class. The .write() method is used.
        @return: nothing
        
        """
        
        regressorNames = self._linearModel.regressorNames()
        if regressorNames[0] == "1":
            outputStream.write("Model: y = a_0")
        else:
            outputStream.write(" + a_0 * %s" % regressorNames[0])
        for n in range(1, self._nParameters):
            if regressorNames[n] == "1":
                outputStream.write(" + a_%d" % n)
            else:
                outputStream.write(" + a_%d * %s" % (n, regressorNames[n]))
        outputStream.write("\n")
        
        coeff = self.regressionCoefficients()
        errors = self.errorBars()
        outputStream.write("Regression coefficients:\n")
        for n in range(self._nParameters):
            outputStream.write("    a_%d = %f +/- %f\n" % (n, coeff[n], errors[n])) 
        
        outputStream.write("Sum of squared residuals: %f\n" % self.sumSqResiduals())
        outputStream.write("Estimated variance of the residuals: %f\n" % self.residualVariance())
        outputStream.write("Coefficient of determination R^2: %f\n" % self.coefficientOfDetermination())
        outputStream.write("F-statistic: %f\n" % self.Fstatistic())
        outputStream.write("AIC: %f\n" % self.AICvalue())
        outputStream.write("BIC: %f\n" % self.BICvalue())



    def __str__(self):
    
        """
        Returns the string written by printing the LinearFit object
         
        """
        
        return   "Linear fit of a dataset of %d observations, " % self._nObservations    \
               + "using a linear model with %d regressors" % self._nParameters











def combinations(items, n):

    """
    Returns all combinations of size n without replacement
    
    @param items: a sequence
    @type items: list, tuple, string, ndarray
    @param n: size of the combinations
    @type n: integer
    @return: generator yielding lists of combinations
    @rtype: generator
    
    Example:
    
    >>> [c for c in combinations("love", 2)]
    [['l', 'o'], ['l', 'v'], ['l', 'e'], ['o', 'v'], ['o', 'e'], ['v', 'e']]
    
    """

    if n==0: 
        yield []
    else:
        for i in xrange(len(items)):
            for cc in combinations(items[i+1:],n-1):
                yield [items[i]]+cc








__all__ = [LinearModel, PolynomialModel, HarmonicModel, LinearFit]
