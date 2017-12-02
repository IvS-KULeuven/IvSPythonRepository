"""
Unit test covering sigproc.fit.py and sigproc.funclib.py

@author: Joris Vos
"""
import copy
import numpy as np
from ivs.sigproc import fit, lmfit, funclib

import unittest
try:
    from mock import patch
    noMock = False
except Exception:
    noMock = True

class FitTestCase(unittest.TestCase):

    def assertArrayEqual(self, l1, l2, msg=None):
        for i, (f1, f2) in enumerate(zip(l1,l2)):
            msg_ = "Array not equal on: %i, %s != %s"%(i, str(f1), str(f2))
            if msg != None: msg_ = msg_ + ", " + msg
            self.assertEqual(f1,f2, msg=msg_)

    def assertArrayAlmostEqual(self, l1,l2,places=None,delta=None, msg=None):
            for i, (f1, f2) in enumerate(zip(l1, l2)):
                msg_ = "Array not equal on: %i, %s != %s"%(i, str(f1), str(f2))
                if msg != None: msg_ = msg_ + ", " + msg
                self.assertAlmostEqual(f1,f2,places=places, delta=delta, msg=msg_)

class TestCase1Function(FitTestCase):

    def setUp(self):
        self.pnames = ['ampl', 'freq', 'phase']
        self.function = lambda p, x: p[0] * np.sin(2*np.pi*(p[1]*x + p[2]))
        self.Dfunction = lambda p, x: p[0] * np.cos(2*np.pi*(p[1]*x + p[2]))
        self.model = fit.Function(function=self.function, par_names=self.pnames,\
                                  jacobian=self.Dfunction)

        self.value = [1.5, 0.12, 3.14]
        self.bounds = [[1.0, 2.0], [-0.20,0.20], [0.0, 6.24]]
        self.vary = [True, True, True]
        self.expression = ['<10', None, None]
        self.model.setup_parameters(values=self.value, bounds=self.bounds, vary=self.vary,\
        exprs=self.expression)

    def testConstructor(self):
        """ sigproc.fit.Function constructor """
        self.assertEqual(self.model.function, self.function)
        self.assertEqual(self.model.jacobian, self.Dfunction)
        self.assertEqual(self.model.par_names, self.pnames)
        self.assertListEqual(self.pnames, self.model.parameters.keys())

    def testEvaluate(self):
        """ sigproc.fit.Function evaluate """
        pars = [1.5, 0.12, 3.14]
        x = np.linspace(0,10,1000)
        y1 = self.function(pars, x)
        y2 = self.model.evaluate(x, pars)
        self.assertTrue(np.all(y1 == y2))

    def testEvaluateJacobian(self):
        """ sigproc.fit.Function evaluate_jacobian """
        pars = [1.5, 0.12, 3.14]
        x = np.linspace(0,10,1000)
        y1 = self.Dfunction(pars, x)
        y2 = self.model.evaluate_jacobian(x, pars)
        self.assertTrue(np.all(y1 == y2))

    def testSetupParameters(self):
        """ sigproc.fit.Function setup_parameters """
        val,err,var,min,max,expr = self.model.get_parameters(full_output=True)
        self.assertArrayEqual(self.value, val)
        self.assertArrayEqual(self.vary, var)
        self.assertArrayEqual(min, np.array(self.bounds)[:,0].tolist())
        self.assertArrayEqual(max, np.array(self.bounds)[:,1].tolist())
        self.assertArrayEqual(expr, self.expression)

    def testUpdateParameters(self):
        """ sigproc.fit.Function  update_parameter """
        self.model.update_parameter('ampl', value=1.75, min=0.0, max=3.5, vary=False)
        val,err,vary,min,max,expr = self.model.get_parameters(parameters=['ampl'], full_output=True)
        self.assertEqual(1.75, val[0])
        self.assertEqual(0.0, min[0])
        self.assertEqual(3.5, max[0])
        self.assertEqual(False, vary[0])

    def testBoundaryAdapting(self):
        """ sigproc.fit.Function setup_parameters boundary adapting """
        pars = [2.0, 0.10, 3.14]
        bounds = [(1.0, 1.5), (0.50, 1.0), (0.0, 6.28)]
        self.model.setup_parameters(values=pars, bounds=bounds)
        val, err, var, min, max, expr = self.model.get_parameters(full_output=True)
        self.assertArrayEqual(min, [1.0, 0.10, 0.0])
        self.assertArrayEqual(max, [2.0, 1.0, 6.28])

class TestCase2Model(unittest.TestCase):

    def setUp(self):
        self.pnames1 = ['p1', 'p2']
        self.function1 = lambda p, x: p[0] + p[1] * x
        self.func1 = fit.Function(function=self.function1, par_names=self.pnames1)

        self.values1 = [2.5, 1.0]
        self.func1.setup_parameters(values=self.values1)

        self.pnames2 = ['p3', 'p4']
        self.function2 = lambda p, x: p[0] * x**2 + p[1] * x**3
        self.func2 = fit.Function(function=self.function2, par_names=self.pnames2)

        self.values2 = [-1.0, 0.5]
        self.func2.setup_parameters(values=self.values2)

        self.mexpr = lambda x: x[0] + x[1]
        self.model = fit.Model(functions=[self.func1, self.func2], expr=self.mexpr)

    def testConstructor(self):
        """ sigproc.fit.Model constructor """
        self.assertListEqual(self.model.functions, [self.func1, self.func2])
        self.assertEqual(self.model.expr, self.mexpr)
        self.assertNotEqual(self.model.par_names, None)

    def testPullParameters(self):
        """ sigproc.fit.Model pull parameters """
        names = self.model.par_names
        self.assertTrue(self.pnames1[0], names[0][0])
        self.assertTrue(self.pnames1[1], names[0][1])
        self.assertTrue(self.pnames2[0], names[1][0])
        self.assertTrue(self.pnames2[1], names[1][1])
        values = self.model.get_parameters()[0].tolist()
        self.assertListEqual(values, np.ravel([self.values1, self.values2]).tolist())

    def testPushParameters(self):
        """ sigproc.fit.Model push parameters """
        self.model.parameters['p1_0'].value = 2.4
        self.assertNotEqual(self.func1.parameters['p1'].value, 2.4)
        self.model.push_parameters()
        self.assertEqual(self.func1.parameters['p1'].value, 2.4)

    def testEvaluate(self):
        """ sigproc.fit.Model evaluate """
        x = np.linspace(0,10,1000)
        y1 = self.model.evaluate(x)
        y2 = self.mexpr([self.function1(self.values1, x), self.function2(self.values2, x)])
        y3 = self.mexpr([self.func1.evaluate(x), self.func2.evaluate(x)])
        self.assertTrue(np.all(y1==y2))
        self.assertTrue(np.all(y1==y3))

    def testUpdateParameters(self):
        """ sigproc.fit.Model update_parameter """
        self.model.update_parameter(parameter='p1_0', value=1.5, min=0.0, max=2.5, expr='<=2.5')
        self.assertEqual(self.model.parameters['p1_0'].value, 1.5)
        self.assertEqual(self.model.parameters['p1_0'].min, 0.0)
        self.assertEqual(self.model.parameters['p1_0'].max, 2.5)
        self.assertEqual(self.model.parameters['p1_0'].expr, '<=2.5')

        with patch.object(self.model, 'push_parameters') as mock:
            self.model.update_parameter(parameter='p1_0', value=2.5)
        mock.assert_called_once()

class TestCase3Minimizer(FitTestCase):
    @classmethod
    def setUpClass(cls):

        parameters = lmfit.Parameters()
        parameters.add('ampl', value=10)
        parameters.add('freq', value=1.)
        parameters.add('phase', value=3.14)
        parameters.add('const', value=-1.)

        model = funclib.sine()

        x = np.linspace(0,10, num=100)
        y = model.evaluate(x, parameters)

        cls.parameters = parameters
        cls.model = model
        cls.x = x
        cls.y = y
        cls.weights = np.ones_like(y)
        cls.errors = np.ones_like(y)

    @unittest.skipIf(noMock, "Mock not installed")
    @patch('ivs.sigproc.lmfit.Minimizer')
    def test1setupResiduals(self, mocked_class):
        """ sigproc.fit.Minimizer _setup_residual_function """

        result = fit.minimize(self.x, self.y, self.model)

        msg = 'attr residuals is not set'
        self.assertTrue(hasattr(result, 'residuals'), msg=msg)

        msg = 'residuals are not correctly calculated'
        residuals = result.residuals(self.parameters, self.x, self.y,
                                     weights=self.weights)
        self.assertTrue(np.all(np.abs(residuals) < 1.0), msg=msg)

    @unittest.skipIf(noMock, "Mock not installed")
    @patch('ivs.sigproc.lmfit.Minimizer')
    def test2setupJacobian(self, mocked_class):
        """ sigproc.fit.Minimizer _setup_jacobian_function """
        result = fit.minimize(self.x, self.y, self.model, err=self.errors)

        msg = 'attr jacobian is not set'
        self.assertTrue(hasattr(result, 'jacobian'), msg=msg)

        msg = 'jacobian should be None in this case'
        self.assertTrue(result.jacobian == None, msg=msg)

    @unittest.skipIf(noMock, "Mock not installed")
    @patch('ivs.sigproc.lmfit.Minimizer')
    def test3perturbdata(self, mocked_class):
        """ sigproc.fit.Minimizer _perturb_input_data """
        result = fit.minimize(self.x, self.y, self.model, errors=self.errors)

        result.x = self.x.reshape(2,50)
        result.y = self.y.reshape(2,50)
        result.errors = result.errors.reshape(2,50)

        np.random.seed(11)
        npoints = 100
        y_ = result._perturb_input_data(100)
        dy = np.abs(np.ravel(y_[0]) - np.ravel(result.y))

        print dy

        msg = 'Perturbed data has wrong shape'
        self.assertEqual(y_.shape, (100,) + result.y.shape, msg=msg)

        msg = 'Perturbed data is NOT perturbed'
        self.assertTrue(np.all(dy > 0.), msg=msg)

        msg = 'Perturbed data is to strongly perturbed'
        self.assertTrue(np.all(dy < 3.), msg=msg)

    @unittest.skipIf(noMock, "Mock not installed")
    @patch('ivs.sigproc.lmfit.Minimizer')
    def test4mcErrorFromParameters(self, mocked_class):
        """ sigproc.fit.Minimizer _mc_error_from_parameters """
        pass


class TestCase4FittingFunction(FitTestCase):
    """
    Integration test testing the fitting of a simple Function
    """

    @classmethod
    def setUpClass(cls):
        #-- setup the function
        pnames = ['ampl', 'freq', 'phase']
        function = lambda p, x: p[0] * np.sin(2*np.pi*(p[1]*x + p[2]))
        Dfunction = lambda p, x: p[0] * np.cos(2*np.pi*(p[1]*x + p[2]))
        fitFunc = fit.Function(function=function, par_names=pnames)
        value = [1.5, 0.12, 3.14]
        bounds = [[1.0, 2.0], [0.0,0.25], [0.0, 6.28]]
        vary = [True, True, True]
        fitFunc.setup_parameters(values=value, bounds=bounds, vary=vary)

        #-- create fake data
        x = np.linspace(0,10,1000)
        np.random.seed(1111)
        y = fitFunc.evaluate(x) + np.random.normal(loc=0.0, scale=0.1, size=1000)

        #-- Change starting parameters
        fitFunc.update_parameter(parameter='ampl', values=1.2)
        fitFunc.update_parameter(parameter='freq', values=0.10)
        fitFunc.update_parameter(parameter='phase', values=3.0)

        cls.x = x
        cls.y = y
        cls.pnames = pnames
        cls.value = value
        cls.fitFunc = fitFunc
        cls.dFunc = Dfunction

    def setUp(self):
        self.model = copy.copy(self.fitFunc)

    def test1minimize(self):
        """ I sigproc.fit.Minimizer Function minimize """
        model = self.model

        result = fit.minimize(self.x, self.y, model)
        values = model.get_parameters()[0]

        msg = 'Fit did not converge to the correct values'
        self.assertAlmostEqual(values[0], self.value[0], places=2, msg=msg)
        self.assertAlmostEqual(values[1], self.value[1], places=2, msg=msg)
        self.assertAlmostEqual(values[2], self.value[2], places=2, msg=msg)

    def test2grid_minimize(self):
        """ I sigproc.fit.Minimizer Function grid_minimize """
        model = self.model
        points = 100
        fitters, newmodels, chisqrs = fit.grid_minimize(self.x, self.y, model,
                                            parameters=self.pnames, points=points,
                                            return_all=True, verbose=False)
        values = newmodels[0].get_parameters()[0]
        pars1 = [mod_.parameters['ampl'] for mod_ in newmodels]

        msg = 'Not correct amount of fitting points'
        self.assertEqual(len(fitters), points, msg=msg)
        self.assertEqual(len(newmodels), points, msg=msg)
        self.assertEqual(len(chisqrs), points, msg=msg)

        msg = 'Parameters not well distributed'
        for i, j in zip(np.random.randint(0,points, size=50), np.random.randint(0,points, size=50)):
            if i != j:
                self.assertNotEqual(pars1[i], pars1[j], msg=msg)

        msg = 'Chi2 is not ordened correctly'
        self.assertTrue(chisqrs[0] == np.min(chisqrs), msg=msg)
        self.assertTrue(np.all(np.diff(chisqrs) >= 0), msg=msg)

        msg = 'Fit did not converge to the correct values'
        self.assertArrayAlmostEqual(values[0:2], self.value[0:2], places=2, msg=msg)

    def test3ci_interval(self):
        """ I sigproc.fit.Minimizer Function calculate_CI """
        model = self.model

        result = fit.minimize(self.x, self.y, model)
        ci = result.calculate_CI(parameters=None, sigmas=[0.674, 0.997])

        expected_ci = "\n\t{'ampl': {0.674: (1.4955, 1.5036), 0.997: (1.4873, 1.5117)}"\
                     +"\n\t'phase': {0.674: (3.1372, 3.1394), 0.997: (3.1349, 3.1417)},"\
                     +"\n\t'freq': {0.674: (0.1202, 0.1206), 0.997: (0.1197, 0.1210)}}"


        msg = 'dict output format is not correct:\nexpected: '+expected_ci+\
              '\nreceived: \n\t'+str(ci)
        self.assertTrue('ampl' in ci, msg=msg)
        self.assertTrue('phase' in ci, msg=msg)
        self.assertTrue('freq' in ci, msg=msg)

        self.assertTrue(0.674 in ci['ampl'] and 0.997 in ci['ampl'], msg=msg)
        self.assertTrue(0.674 in ci['phase'] and 0.997 in ci['phase'], msg=msg)
        self.assertTrue(0.674 in ci['freq'] and 0.997 in ci['freq'], msg=msg)

        msg = 'ci results are not correct, expected: \n' + expected_ci
        self.assertArrayAlmostEqual(ci['ampl'][0.674], [1.4955, 1.5036], places=2, msg=msg)
        self.assertArrayAlmostEqual(ci['ampl'][0.997], [1.4873, 1.5117], places=2, msg=msg)
        self.assertArrayAlmostEqual(ci['phase'][0.674], [3.1372, 3.1394], places=2, msg=msg)
        self.assertArrayAlmostEqual(ci['phase'][0.997], [3.1349, 3.1417], places=2, msg=msg)
        self.assertArrayAlmostEqual(ci['freq'][0.674], [0.1202, 0.1206], places=2, msg=msg)
        self.assertArrayAlmostEqual(ci['freq'][0.997], [0.1197, 0.1210], places=2, msg=msg)

        msg = 'ci intervals are not stored in the parameter object'
        for name, param in model.parameters.items():
            self.assertTrue(0.674 in param.cierr, msg=msg)
            self.assertTrue(0.997 in param.cierr, msg=msg)

        msg = 'Reruning calculate_CI fails and raises an exception!'
        try:
            result.calculate_CI(parameters=None, sigma=0.85, maxiter=200, short_output=True)
            result.calculate_CI(parameters=None, sigma=0.95, maxiter=200, short_output=True)
        except Exception:
            self.fail(msg)

    def test4ci2d_interval(self):
        """ I sigproc.fit.Minimizer Function calculate_CI_2D """
        model = self.model

        result = fit.minimize(self.x, self.y, model)
        x, y, ci = result.calculate_CI_2D(xpar='ampl', ypar='freq', res=10, limits=None,
                                          ctype='prob')
        values, errors = model.get_parameters()

        msg = 'Auto limits are not value +- 5 * stderr'
        vmin = values - 5*errors
        vmax = values + 5*errors
        self.assertAlmostEqual(min(x), vmin[0], places=4, msg=msg)
        self.assertAlmostEqual(max(x), vmax[0], places=4, msg=msg)
        self.assertAlmostEqual(min(y), vmin[1], places=4, msg=msg)
        self.assertAlmostEqual(max(y), vmax[1], places=4, msg=msg)

        msg = 'Grid values are not correctly distributed'
        exp1 = [1.4790, 1.4835, 1.4881, 1.4927, 1.4972,
                1.5018, 1.5064, 1.5109, 1.5155, 1.5200]
        exp2 = [0.11937, 0.11960, 0.11983, 0.12006, 0.12028,
                0.12051, 0.12074, 0.12097, 0.12120, 0.12143]
        self.assertArrayAlmostEqual(x, exp1, places=3, msg=msg)
        self.assertArrayAlmostEqual(y, exp2, places=4, msg=msg)

        msg = 'resolution of the grid is not correct'
        self.assertEqual(len(ci), 10, msg=msg)
        self.assertEqual(len(ci[0]), 10, msg=msg)

        msg = 'CI values are not correct'
        exp1 = [99.9996, 99.9572, 98.2854, 79.4333, 27.5112,
                25.6265, 77.7912, 98.0525, 99.9489, 99.9995]
        exp2 = [99.9996, 99.9565, 98.2751, 79.4030, 27.5112,
                25.6436, 77.8194, 98.0628, 99.9496, 99.9996]
        self.assertArrayAlmostEqual(ci[4], exp1, places=2, msg=msg)
        self.assertArrayAlmostEqual(ci[:,4], exp2, places=2, msg=msg)

    def test5mc_error(self):
        """ I sigproc.fit.Minimizer Function calculate_MC_error """
        result = fit.minimize(self.x, self.y, self.model)

        np.random.seed(11)
        mcerrors = result.calculate_MC_error(errors=0.5, points=100, short_output=True,
                                             verbose=False)

        msg = 'Short output returned wrong type (array expected)'
        self.assertEqual(np.shape(mcerrors), (3,), msg=msg)

        msg = 'Calculated MC errors are wrong'
        self.assertArrayAlmostEqual(mcerrors, [0.02156, 0.00108, 0.00607], places=3, msg=msg)

        np.random.seed(11)
        mcerrors = result.calculate_MC_error(errors=0.5, points=50, short_output=False,
                                             verbose=False)

        msg = 'Long output returned wrong type (dict expected)'
        self.assertEqual(type(mcerrors), dict, msg=msg)

        msg = 'Not all parameters are included in output'
        self.assertTrue('ampl' in mcerrors, msg=msg)
        self.assertTrue('freq' in mcerrors, msg=msg)
        self.assertTrue('phase' in mcerrors, msg=msg)

        msg = 'MC error is not stored in the parameters object'
        params = self.model.parameters
        self.assertEqual(params['ampl'].mcerr, mcerrors['ampl'], msg=msg)
        self.assertEqual(params['freq'].mcerr, mcerrors['freq'], msg=msg)
        self.assertEqual(params['phase'].mcerr, mcerrors['phase'], msg=msg)

class TestCase5IntegrationJacobian(FitTestCase):
    """
    Integration test testing the fitting of a Function with Jacobian
    """
    @classmethod
    def setUpClass(cls):
        "===Exponential_decay_Jacobian==="

        pnames = ['p1', 'p2', 'p3']

        def func(var, x):
            return var[0]*np.exp(-var[1]*x)+var[2]

        def dfunc(var,x):
            v = np.exp(-var[1]*x)
            return np.array([-v,+var[0]*x*v,-np.ones(len(x))]).T

        cls.function = func
        cls.jacobian = dfunc

        model = fit.Function(function=func, par_names=pnames, jacobian=dfunc)
        pars = [10, 10, 10]
        bounds = [(0, 20), (0, 20), (0, 20)]
        vary = [True, True, True]
        model.setup_parameters(values=pars, bounds=bounds, vary=vary)
        cls._model = model

        #-- Creating data to fit
        xs = np.linspace(0,4,50)
        ys = func([2.5,1.3,0.5],xs)
        np.random.seed(3333)
        err = 0.1*np.random.normal(size=len(xs))
        yn = ys + err

        cls.x = xs
        cls.y = yn
        cls.err = np.abs(err)

    def setUp(self):
        self.model = copy.deepcopy(self._model)

    def test1minimize(self):
        """ I sigproc.fit.Minimizer Jacobian minimize """

        self.model.setup_parameters(value=[3.0, 1.0, 1.0])
        result = fit.minimize(self.x, self.y, self.model)
        names = self.model.parameters.name
        val = self.model.parameters.value
        err = self.model.parameters.stderr
        valr = [2.4272, 1.2912, 0.4928]
        errr = [0.0665, 0.0772, 0.0297]

        print self.model.param2str(accuracy=4, output='result')

        msg = 'Resulting value for {:s} is not correct within {:.0f}\%'
        pc = 0.05 # 5% difference is allowed
        for i in range(3):
            self.assertAlmostEqual(val[i], valr[i], delta = pc*valr[i],
                                   msg=msg.format(names[i], pc*100))

        msg = 'Resulting std-error for {:s} is not correct within {:.0f}\%'
        pc = 0.1 # 10% difference is allowed
        for i in range(3):
            self.assertAlmostEqual(err[i], errr[i], delta = pc*errr[i],
                                   msg=msg.format(names[i], pc*100))

    def test2grid_minimize(self):
        """ I sigproc.fit.Minimizer Jacobian grid_minimize """

        result = fit.grid_minimize(self.x, self.y, self.model, points=100)
        names = self.model.parameters.name
        val = self.model.parameters.value
        err = self.model.parameters.stderr
        valr = [2.4272, 1.2912, 0.4928]
        errr = [0.4345, 0.3796, 0.0921]

        print self.model.param2str(accuracy=4, output='result')

        msg = 'Resulting value for {:s} is not correct within {:.0f}\%'
        pc = 0.05 # 5% difference is allowed
        for i in range(3):
            self.assertAlmostEqual(val[i], valr[i], delta = pc*valr[i],
                                   msg=msg.format(names[i], pc*100))

        msg = 'Resulting std-error for {:s} is not correct within {:.0f}\%'
        pc = 0.1 # 10% difference is allowed
        for i in range(3):
            self.assertAlmostEqual(err[i], errr[i], delta = pc*errr[i],
                                   msg=msg.format(names[i], pc*100))

    def test3ci_interval(self):
        """ I sigproc.fit.Minimizer Jacobian calculate_CI """
        self.model.setup_parameters(value=[3.0, 1.0, 1.0])

        result = fit.minimize(self.x, self.y, self.model)
        ci = result.calculate_CI(parameters=None, sigmas=[0.994])
        cir = [[2.2378,2.6210],[1.0807,1.5276],[0.4015,0.5740]]
        val = self.model.parameters.value
        names = self.model.parameters.name

        msg = 'The 99.4\% CI interval for {:s} is not correct within {:.0f}\%'
        pc = 0.05
        for i,p in enumerate(names):
            self.assertArrayAlmostEqual(ci[p][0.994], cir[i], delta = pc*val[i],
                                   msg=msg.format(p, pc*100))

    def test4ci2d_interval(self):
        """ I sigproc.fit.Minimizer Jacobian calculate_CI_2D """
        self.model.setup_parameters(value=[3.0, 1.0, 1.0])

        result = fit.minimize(self.x, self.y, self.model)
        x, y, ci = result.calculate_CI_2D(xpar='p1', ypar='p2', res=10, limits=None,
                                          ctype='prob')

        cir = [99.9977, 99.9019, 97.5451, 74.5012, 20.5787,
               38.3632, 87.5572, 99.1614, 99.9734, 99.9994]

        msg = 'The 2D CI values are not correct within 2\%'
        pc = 0.05
        for i in range(3):
            self.assertArrayAlmostEqual(ci[4], cir, delta = 2, msg = msg)

    def test5mc_error(self):
        """ I sigproc.fit.Minimizer Jacobian calculate_MC_error """
        self.model.setup_parameters(value=[3.0, 1.0, 1.0])
        result = fit.minimize(self.x, self.y, self.model)
        result.calculate_MC_error(errors=self.err, points=100, short_output=True,
                                  verbose=False)

        names = self.model.parameters.name
        err = self.model.parameters.mcerr
        errr = [0.0795, 0.0688, 0.0281]

        msg = 'Resulting MC-error for {:s} is not correct within {:.0f}\%'
        pc = 0.05 # 10% difference is allowed
        for i in range(3):
            self.assertAlmostEqual(err[i], errr[i], delta = pc*errr[i],
                                   msg=msg.format(names[i], pc*100))



class TestCase6FittingModel(FitTestCase):
    """
    Integration test testing the fitting of a simple Model
    """

    @classmethod
    def setUpClass(cls):
        gauss1 = funclib.gauss()
        pars = [-0.75,1.0,0.1,1]
        gauss1.setup_parameters(values=pars)
        gauss2 = funclib.gauss()
        pars = [0.22,1.0,0.01,0.0]
        vary = [True, True, True, False]
        gauss2.setup_parameters(values=pars, vary=vary)

        mymodel = fit.Model(functions=[gauss1, gauss2])

        np.random.seed(1111)
        x = np.linspace(0.5,1.5, num=1000)
        y = mymodel.evaluate(x)
        noise = np.random.normal(0.0, 0.015, size=len(x))
        y = y+noise

        pars = [-0.70,1.0,0.11,0.95]
        gauss1.setup_parameters(values=pars)
        pars = [0.27,1.0,0.005,0.0]
        vary = [True, True, True, False]
        gauss2.setup_parameters(values=pars, vary=vary)

        cls.x = x
        cls.y = y
        cls.value1 = [-0.75,1.0,0.1,1]
        cls.value2 = [0.22,1.0,0.01,0.0]
        cls.fitModel = mymodel

    def testMinimize(self):
        """ I sigproc.fit.Minimizer Model minimize """
        model = self.fitModel

        result = fit.minimize(self.x, self.y, model)
        functions = model.functions
        values1 = functions[0].get_parameters()[0]
        values2 = functions[1].get_parameters()[0]

        self.assertArrayAlmostEqual(values1, self.value1, places=2)
        self.assertArrayAlmostEqual(values2, self.value2, places=2)

class TestCase7Funclib(unittest.TestCase):

    def testFunclibFunctions(self):
        """ sigproc.funclib: all functions """
        def check_function(fname, param,  x, expected, dec=5, **kwargs):
            y = funclib.evaluate(fname,np.array(x),param, **kwargs)
            y = np.round(y, dec)
            self.assertTrue(np.all(y == expected), msg=fname+" Failed")

        #-- blackbody
        check_function('blackbody', [5798, 1e-3], [3000, 4120, 10734], [3141.80541, 6106.61968, 2304.17307])

        #-- rayleigh_jeans
        check_function('rayleigh_jeans', [5798, 1e-5], [3000, 4120, 10734], [14853.11541, 4175.55021, 90.62671])

        #-- wien
        check_function('wien', [5798, 1e-3], [3000, 4120, 10734], [3141.00219, 6091.82801, 2075.87277])

        #-- box_transit
        check_function('box_transit', [2.,0.4,0.1,0.3,0.5], [0.74, 0.75, 2.74, 2.75], [1.5, 2.0, 2.0, 1.5])

        #-- polynomial
        check_function('polynomial', [1,2.5,1.2], [0.27, 1.34, 2.01], [1.9479, 6.3456, 10.2651], d=2)

        #-- soft_parabola
        check_function('soft_parabola', [1.,0,1.,0.25], [-0.87, 0, 1.0], [0.83796, 1.0, 0.0])

        #-- gauss
        check_function('gauss', [5,1.,2.,0.5], [-5, 1.0, 4.23], [0.55554, 5.5, 1.85707])

        #-- voigt
        check_function('voigt', [20.,1.,1.5,3.,0.5], [-5, 1.0, 4.23], [0.97541, 2.28835, 1.5696])

        #-- lorentz
        check_function('lorentz', [5,1.,2.,0.5], [-5, 1.0, 4.23], [0.625, 1.75, 0.84643])

        #-- sine
        check_function('sine', [1.,2.,0,0], [0.24, 2.0, 7.32], [0.12533, 0.0, -0.77051])

        #-- sine_linfreqshift
        check_function('sine_linfreqshift', [1.,0.5,0,0,.5], [0.24, 2.0, 7.32], [0.74761, 0.0, 0.34228])

        #-- sine_expfreqshift
        check_function('sine_expfreqshift', [1.,0.5,0,0,1.2], [0.24, 2.0, 7.32], [0.69665, 0.96315, -0.88957])

        #-- sine_orbit
        check_function('sine_orbit', [1.,2.,0,0,0.1,10.,0.1], [0.24, 2.0, 7.32], [0.01663, 0.63673, -0.15785])

        #-- kepler_orbit
        check_function('kepler_orbit', [2.5,0.,0.5,0,3,1.], [0.24, 2.0, 7.32], [2.55123, 0.63411, 3.34141])

        #-- power_law
        check_function('power_law', [2.,3.,1.5,0.0,0.5], [0.24, 2.0, 7.32], [1.74151, 0.62741, 0.51925])

#if __name__ == '__main__':
    #unittest.main()


#suite = [unittest.TestLoader().loadTestsFromTestCase(LMfitTestCase)]
#suite += unittest.TestLoader().loadTestsFromTestCase(FunctionTestCase)
#allTests = unittest.TestSuite(suite)
#unittest.TextTestRunner(verbosity=2).run(allTests)




















