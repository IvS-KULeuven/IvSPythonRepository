"""
Unit test covering sigproc.lmfit

@author: Joris Vos
"""
import copy
import numpy as np
from ivs.sigproc.lmfit import Minimizer, ConfidenceInterval
from ivs.sigproc.lmfit import Parameter
from ivs.sigproc import lmfit

import unittest
try:
    import mock
    from mock import MagicMock
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

class TestCase1Parameter(FitTestCase):

    def setUp(self):
        self.attrs1 = dict(name='p1', value=15.21, vary=True, min=10., max=20., expr=None)
        self.attrs2 = dict(name='p2', value=-0.22, vary=True, min=None, max=None, expr=None)
        self.param1 = Parameter(**self.attrs1)
        self.param2 = Parameter(**self.attrs2)

        self.param1.stderr = 0.1521
        self.param1.mcerr = 0.7605
        self.ci = {'0.95':(14.10, 16.32), '0.997':(13.52, 17.67)}
        self.param1.cierr = copy.copy(self.ci)

    def test1NewAttrs(self):
        """ sigproc.lmfit.Parameter: new error attributes """
        par = self.param1
        self.assertTrue(hasattr(par, 'stderr'), msg='Parameter does not have var stderr')
        self.assertTrue(hasattr(par, 'correl'), msg='Parameter does not have var correl')
        self.assertTrue(hasattr(par, 'cierr'), msg='Parameter does not have var cierr')
        self.assertTrue(hasattr(par, 'mcerr'), msg='Parameter does not have var mcerr')

    def test2CanKick(self):
        """ sigproc.lmfit.Parameter: can_kick() """
        par1 = self.param1
        par2 = self.param2
        self.assertTrue(par1.can_kick(), msg='Parameter should be able to kick')
        self.assertFalse(par2.can_kick(), msg='Parameter should NOT be able to kick')

    def test3Kick(self):
        """ sigproc.lmfit.Parameter: kick() """
        par = self.param1
        par.kick()

        msg = 'Parameter value needs to change after Kick'
        self.assertNotEqual(par.value, self.attrs1['value'], msg=msg)
        self.assertEqual(par.value, par.user_value, msg=msg)

        msg = 'Kicked Parameter value should be between min and max'
        self.assertTrue(par.value >= self.attrs1['min'] and par.value <= self.attrs1['max'],
                        msg=msg)

    def test4GetAttr(self):
        """ sigproc.lmfit.Parameter: __getattr__() """
        par = self.param1

        msg = 'CI values should be reachable directly by fx. ci95'
        self.assertTrue(hasattr(par, 'ci95'), msg=msg)
        self.assertTrue(hasattr(par, 'ci997'), msg=msg)
        self.assertTrue(hasattr(par, 'ci665'), msg=msg)

        msg = 'Returned CI values are not correct'
        self.assertEqual(par.ci95, self.ci['0.95'], msg=msg)
        self.assertEqual(par.ci997, self.ci['0.997'], msg=msg)
        self.assertEqual(par.ci665, (np.nan, np.nan), msg=msg)

        msg = 'Percentual error should be reachable directly by fx. stderrpc'
        self.assertTrue(hasattr(par, 'stderrpc'), msg=msg)
        self.assertTrue(hasattr(par, 'mcerrpc'), msg=msg)

        msg = 'Returned percentual errors are not correct'
        self.assertEqual(par.stderrpc, 1., msg=msg)
        self.assertEqual(par.mcerrpc, 5., msg=msg)

class TestCase2Parameters(FitTestCase):

    def setUp(self):
        self.attrs1 = dict(value=15.21, vary=True, min=10., max=20., expr=None)
        self.attrs2 = dict(value=-0.22, vary=True, min=None, max=None, expr=None)
        self.attrs3 = dict(value=None, vary=True, min=None, max=None, expr='p2/2.0')
        self.attrs4 = dict(value=1000., vary=True, min=100, max=1500, expr=None)

        self.params = lmfit.Parameters()
        self.params.add('p1', **self.attrs1)
        self.params.add('p2', **self.attrs2)
        self.params.add('p3', **self.attrs3)
        self.params.add('p4', **self.attrs4)

    def test1CanKick(self):
        """ sigproc.lmfit.Parameters: can_kick() """

        msg='can_kick() does not return the correct parameters (p1, p4)'
        kick_pars = self.params.can_kick()
        self.assertEqual(kick_pars, ['p1', 'p4'], msg=msg)

        msg='can_kick(pnames) does not return the correct parameters (p4)'
        kick_pars = self.params.can_kick(pnames=['p2', 'p3', 'p4'])
        self.assertEqual(kick_pars, ['p4'], msg=msg)

    @unittest.skipIf(noMock, "Mock not installed")
    def test2Kick(self):
        """ sigproc.lmfit.Parameters: kick() """

        self.params['p1'].kick = MagicMock()
        self.params['p2'].kick = MagicMock()
        self.params['p3'].kick = MagicMock()
        self.params['p4'].kick = MagicMock()

        self.params.kick(pnames=['p1'])
        msg='kick(pnames) did not kick the requested parameter (p1)'
        self.assertTrue(self.params['p1'].kick.called, msg=msg)
        msg='kick(pnames) should NOT kick other parameters (p2, p3, p4)'
        self.assertFalse(self.params['p2'].kick.called, msg=msg)
        self.assertFalse(self.params['p3'].kick.called, msg=msg)
        self.assertFalse(self.params['p4'].kick.called, msg=msg)

        self.params['p1'].kick = MagicMock()
        self.params['p2'].kick = MagicMock()
        self.params['p3'].kick = MagicMock()
        self.params['p4'].kick = MagicMock()

        self.params.kick()
        msg='kick() did not kick all possible parameters (p1, p4)'
        self.assertTrue(self.params['p1'].kick.called, msg=msg)
        self.assertTrue(self.params['p4'].kick.called, msg=msg)
        msg='kick() should NOT kick other parameters (p2, p3)'
        self.assertFalse(self.params['p2'].kick.called, msg=msg)
        self.assertFalse(self.params['p3'].kick.called, msg=msg)

    def test3test4GetAttr(self):
        """ sigproc.lmfit.Parameters: __getattr__() """

        pars = self.params

        msg = 'Parameters does not have attribute {:s}'
        for name in ['name', 'value', 'min', 'max', 'vary', 'expr',
                     'stderr', 'mcerr', 'cierr']:
            self.assertTrue(hasattr(pars, name), msg=msg.format(name))

        msg = 'value attr in parameters does NOT return the correct values'
        self.assertArrayEqual(pars.value, [15.21, -0.22, None, 1000], msg=msg)

        msg = 'value min in parameters does NOT return the correct values'
        self.assertArrayEqual(pars.min, [10., -np.inf, -np.inf, 100], msg=msg)

        msg = 'value expr in parameters does NOT return the correct values'
        self.assertArrayEqual(pars.expr, [None, None, 'p2/2.0', None], msg=msg)

class TestCase3ConfidenceIterval(FitTestCase):

    @unittest.skipIf(noMock, "Mock not installed")
    @mock.patch.object(ConfidenceInterval, '__init__')
    def test1Format(self, mocked_init):
        """ sigproc.lmfit.confidence.conf_interval(): output format """
        mocked_init.return_value = None

        ci = ConfidenceInterval()
        ci.p_names = ['p1', 'p2']
        ci.sigmas = [0.674, 0.95, 0.997]
        ci.trace = False
        ci.calc_ci = MagicMock(return_value=[('0.674', 10.0), ('0.95',-0.022), ('0.997', 56.)])
        ciOut =  ci.calc_all_ci()

        expected_output = "{'p1': {0.95: (.., ..), 0.674: (.., ..), 0.997: (.., ..)}" + \
                          ", 'p2': {0.95: (.., ..), 0.674: (.., ..), 0.997: (.., ..)}}"

        msg = 'Both parameters (p1, p2) should be in the output dict.' + expected_output
        self.assertTrue('p1' in ciOut, msg=msg)
        self.assertTrue('p2' in ciOut, msg=msg)

        msg = 'All sigma values should be in the output dictionary.' + expected_output
        self.assertTrue(0.674 in ciOut['p1'], msg=msg)
        self.assertTrue(0.95 in ciOut['p1'], msg=msg)
        self.assertTrue(0.997 in ciOut['p1'], msg=msg)
        self.assertTrue(0.674 in ciOut['p2'], msg=msg)
        self.assertTrue(0.95 in ciOut['p2'], msg=msg)
        self.assertTrue(0.997 in ciOut['p2'], msg=msg)

class TestCase4Minimizer(FitTestCase):

    @unittest.skipIf(noMock, "Mock not installed")
    @mock.patch.object(Minimizer, '__init__')
    def test1start_minimizer(self, mocked_init):
        """ sigproc.lmfit.minimize.Minimizer(): start_minimize() """
        mocked_init.return_value = None

        mini = Minimizer()
        mini.anneal = MagicMock(return_value=None)
        mini.lbfgsb = MagicMock(return_value=None)
        mini.fmin = MagicMock(return_value=None)
        mini.leastsq = MagicMock(return_value=None)

        msg = "start_minimize('leastsq') should only start minimizer.leastsq()"
        mini.start_minimize('leastsq')
        self.assertTrue(mini.leastsq.called, msg=msg)
        self.assertFalse(mini.anneal.called, msg=msg)
        self.assertFalse(mini.lbfgsb.called, msg=msg)
        self.assertFalse(mini.fmin.called, msg=msg)
        mini.leastsq = MagicMock(return_value=None)

        msg = "start_minimize('leastsq') should only start minimizer.leastsq()"
        mini.start_minimize('nelder')
        self.assertFalse(mini.leastsq.called, msg=msg)
        self.assertFalse(mini.anneal.called, msg=msg)
        self.assertFalse(mini.lbfgsb.called, msg=msg)
        self.assertTrue(mini.fmin.called, msg=msg)

class TestCase5Integration(FitTestCase):

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
        cls.p_true = p_true

        n = 2500
        xmin = 0.
        xmax = 250.0
        np.random.seed(11)
        noise = np.random.normal(scale=0.7215, size=n)
        x     = np.linspace(xmin, xmax, n)
        data  = residual(p_true, x) + noise

        # setup reasonable starting values for the fit
        p_fit = lmfit.Parameters()
        p_fit.add('amp', value=13.0)
        p_fit.add('period', value=2)
        p_fit.add('shift', value=0.0)
        p_fit.add('decay', value=0.02)

        result = lmfit.minimize(residual, p_fit, args=(x,), kws={'data':data})
        cls.result = result
        cls.p_fit = copy.deepcopy(p_fit)
        ci = lmfit.conf_interval(result)
        cls.ci = ci

    def test1minimizer_value(self):
        """ I sigproc.lmfit minimize: value"""

        p_fit = self.p_fit
        names = p_fit.name
        values = p_fit.value

        p_true = self.p_true
        real = p_true.value

        msg = 'Resulting value for {:s} is not correct within {:.0f}\%'
        pc = 0.05 # 5% difference is allowed
        for i in range(4):
            self.assertAlmostEqual(values[i], real[i], delta = pc*real[i],
                                   msg=msg.format(names[i], pc*100))


    def test2minimizer_stderr(self):
        """ I sigproc.lmfit minimize: stderror"""

        p_fit = self.p_fit
        names = p_fit.name
        err = p_fit.stderr

        msg = 'Resulting std-error for {:s} is not correct within {:.0f}\%'
        pc = 0.1 # 10% difference is allowed
        real = [0.0501, 0.0027, 0.0049, 4.1408e-05]
        for i in range(4):
            self.assertAlmostEqual(err[i], real[i], delta = pc*real[i],
                                   msg=msg.format(names[i], pc*100))

    def test3minimizer_correl(self):
        """ I sigproc.lmfit minimize: correlation """

        p_fit = self.p_fit
        corr = p_fit.correl

        msg = 'Correlations for {:s} is not correct within 10\%'
        self.assertAlmostEqual(corr[1]['shift'], 0.800155, delta = 0.080,
                               msg = msg.format('period-shift'))
        self.assertAlmostEqual(corr[3]['amp'], 0.57584, delta = 0.057,
                               msg = msg.format('amp-decay'))

        msg = 'Correlation for {:s} is to high (> 0.10)'
        self.assertTrue(corr[0]['period'] < 0.10, msg = msg.format('amp-period'))
        self.assertTrue(corr[1]['decay'] < 0.10, msg = msg.format('shift-decay'))
        self.assertTrue(corr[2]['decay'] < 0.10, msg = msg.format('period-decay'))

    def test4confidenceIntervals(self):
        """ I sigproc.lmfit conf_interval: values """
        ci = self.ci

        msg = "{:.0f}\% CI values for '{:s}' not correct"
        self.assertArrayAlmostEqual(ci['period'][0.674], (5.3240, 5.3294), places=2,
                                    msg = msg.format(67.4, 'period'))
        self.assertArrayAlmostEqual(ci['amp'][0.95], (13.8623, 14.0592), places=2,
                                    msg = msg.format(95, 'amp'))
        self.assertArrayAlmostEqual(ci['shift'][0.997], (0.1059, 0.1351), places=2,
                                    msg = msg.format(99.7, 'shift'))


