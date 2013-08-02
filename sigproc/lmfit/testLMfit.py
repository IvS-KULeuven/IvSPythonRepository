"""
Unit test covering sigproc.lmfit
checking that the lmfit package is working and compatible with sigproc.fit

@author: Joris Vos
"""
import copy
import unittest
import numpy as np
from ivs.sigproc import lmfit

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

class TestCase1LMFit(FitTestCase):
    
    @classmethod  
    def setUpClass(cls):
        # the fitting model (residuals)
        def residual(pars, x, data=None):
            amp = pars['amp'].value
            per = pars['period'].value
            shift = pars['shift'].value
            decay = pars['decay'].value

            if abs(shift) > np.pi/2:
                shift = shift - np.sign(shift)*np.pi
            model = amp*np.sin(shift + x/per) * np.exp(-x*x*decay*decay)
            if data is None:
                return model
            return (model - data)
        
        # create sample data
        p_true = lmfit.Parameters()
        p_true.add('amp', value=14.0)
        p_true.add('period', value=5.33)
        p_true.add('shift', value=0.123)
        p_true.add('decay', value=0.010)
        
        n = 2500
        xmin = 0.
        xmax = 250.0
        np.random.seed(11)
        noise = np.random.normal(scale=0.7215, size=n)
        x     = np.linspace(xmin, xmax, n)
        data  = residual(p_true, x) + noise
        
        # setup reasonable starting values for the fit
        fit_params = lmfit.Parameters()
        fit_params.add('amp', value=13.0)
        fit_params.add('period', value=2)
        fit_params.add('shift', value=0.0)
        fit_params.add('decay', value=0.02)
        
        out = lmfit.minimize(residual, fit_params, args=(x,), kws={'data':data})
        cls.fit_params = copy.deepcopy(fit_params)
        ci = lmfit.conf_interval(out)
        cls.ci = ci
        
    
    def setUp(self):
        self.pars = lmfit.Parameters()
        self.pars.add('test1', value=15, min=10, max=20)
        self.pars.add('test2', value=2.39, min=-4.5, max=4.5)
        
        def f(var,xs):
            return var[0]*np.exp(-var[1]*xs)+var[2]     
        def res(par, x, data=None):
            p1 = par['p1'].value
            p2 = par['p2'].value
            p3 = par['p3'].value
            
            if data == None:
                return f([p1,p2,p3], x)
            return data - f([p1,p2,p3], x)
            
        params = lmfit.Parameters()
        params.add('p1', value=10)
        params.add('p2', value=10)
        params.add('p3', value=10)

        xs = np.linspace(0,4,50)
        np.random.seed(3333)
        ys = f([2.5,1.3,0.5],xs) + 0.1*np.random.normal(size=len(xs))
        
        result = lmfit.minimize(res, params, args=(xs, ys))
        #self.ci = lmfit.conf_interval(result, p_names=['p1','p2','p3'], sigmas=[0.95, 0.99])
    
    def test1LMfitMinimize(self):
        """ sigproc.lmfit minimize: value and stderror """
        
        pars = self.fit_params
        values = [pars[p].value for p in pars.keys()]
        err = [pars[p].stderr for p in pars.keys()]
        corr = [pars[p].correl for p in pars.keys()]
        
        self.assertArrayAlmostEqual(values, [13.9606, 5.3267, 0.1205, 0.0099], places=2, 
                                    msg='Fit values are not correct')
        self.assertArrayAlmostEqual(err, [0.0501, 0.0027, 0.0049, 4.1408e-05], places=3,
                                    msg='Fit errors are not correct')
    
    def test2LMfitMinimize(self):
        """ sigproc.lmfit minimize: correlation """
        pars = self.fit_params
        corr = {}
        for p in pars.keys():
            corr[p] = pars[p].correl
        
        self.assertAlmostEqual(corr['period']['shift'], 0.800, places=2, 
                               msg='Correlation period-shift is wrong')
        self.assertAlmostEqual(corr['amp']['decay'], 0.576, places=2,
                               msg='Correlation amp-decay is wrong')
        self.assertTrue(corr['amp']['period'] < 0.10, msg='Correlation amp-period is to high')
        self.assertTrue(corr['shift']['decay'] < 0.10, msg='Correlation shift-decay is to high')
        self.assertTrue(corr['period']['decay'] < 0.10, msg='Correlation period-decay is to high')
    
    def test3ParameterKicking(self):
        """ sigproc.lmfit Parameter Kicking """
        self.assertTrue(hasattr(self.pars, 'kick'))
        self.assertTrue(hasattr(self.pars['test1'], 'kick'))
        self.pars.kick()
        self.assertGreaterEqual(self.pars['test1'].value, self.pars['test1'].min)
        self.assertGreaterEqual(self.pars['test2'].value, self.pars['test2'].min)
        self.assertLessEqual(self.pars['test1'].value, self.pars['test1'].max)
        self.assertLessEqual(self.pars['test2'].value, self.pars['test2'].max)
        
    def test4ConfidenceIntervallOutputFormat(self):
        """ sigproc.lmfit conf_interval: format """
        self.assertTrue('amp' in self.ci, msg='parameters not correct in ci dictionary')
        self.assertTrue('decay' in self.ci, msg='parameters not correct in ci dictionary')
        self.assertTrue('period' in self.ci, msg='parameters not correct in ci dictionary')
        self.assertTrue('shift' in self.ci, msg='parameters not correct in ci dictionary')
        
        self.assertTrue(0.95 in self.ci['period'], msg='sigmas not correct in ci dictionary')
        self.assertTrue(0.674 in self.ci['decay'], msg='sigmas not correct in ci dictionary')
        self.assertTrue(0.997 in self.ci['shift'], msg='sigmas not correct in ci dictionary')
        
        self.assertEqual(len(self.ci['amp'][0.95]), 2, msg='number of ci values not correct')
        
    def test5ConfidenceIntervallOutputValues(self):
        """ sigproc.lmfit conf_interval: values """
        ci = self.ci
        
        self.assertArrayAlmostEqual(ci['period'][0.674], (5.3240, 5.3294), places=2,
                                    msg='CI values not correct')
        self.assertArrayAlmostEqual(ci['amp'][0.95], (13.8623, 14.0592), places=2,
                                    msg='CI values not correct')
        self.assertArrayAlmostEqual(ci['shift'][0.997], (0.1059, 0.1351), places=2,
                                    msg='CI values not correct')
        
        
        
        