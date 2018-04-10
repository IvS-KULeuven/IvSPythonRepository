"""
Unit test covering sed.fit.py

@author: Joris Vos
"""
import numpy as np
from numpy import inf, array
from ivs import sigproc
from ivs.sed import fit, model, builder, filters
from ivs.units import constants
from ivs.catalogs import sesame
from matplotlib import mlab

import unittest
try:
    import mock
    from mock import patch
    noMock = False
except Exception:
    noMock = True

# Set at True to skip integration tests.
noIntegration = False

class SEDTestCase(unittest.TestCase):
    """Add some extra usefull assertion methods to the testcase class"""

    def create_patch(self, name, method, **kwargs):
        patcher = patch.object(name, method, **kwargs)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing

    def assertArrayEqual(self, l1, l2, msg=None):
        msg_ = "Array is not equal to expected array: %s != %s"%(l1,l2)
        if msg != None: msg_ = msg
        self.assertEqual(len(l1), len(l2), msg=msg_)
        self.assertEqual(l1.dtype, l2.dtype, msg=msg_)
        for f1, f2 in zip(l1,l2):
            self.assertEqual(f1,f2, msg=msg_)

    def assertArrayAlmostEqual(self, l1,l2,places=None,delta=None, msg=None):
            msg_ = "assertArrayAlmostEqual Failed: %s != %s"%(str(l1),str(l2))
            if msg != None: msg_ = msg_ + msg
            for f1, f2 in zip(l1, l2):
                self.assertAlmostEqual(f1,f2,places=places, delta=delta, msg=msg_)

    def assertAllBetween(self, l, low, high, msg=None):
        self.assertGreaterEqual(min(l), low, msg=msg)
        self.assertLessEqual(max(l), high, msg=msg)

    def assertAllBetweenDiff(self, l, low, high, places=4, msg=None):
        self.assertGreaterEqual(min(l), low, msg=msg)
        self.assertLessEqual(max(l), high, msg=msg)
        lr = np.round(l, decimals=places)
        self.assertTrue(len(np.unique(lr)) > 1, msg="All elements are equal: %s"%(l))

    def assertInList(self, el, lst, msg=None):
        msg_ = "The given element %s is not in list %s"%(el, lst)
        if msg != None: msg_=msg
        res = []
        for l in lst:
            try:
                res.append(all(el == l))
            except Exception:
                res.append(el == l)
        self.assertTrue(any(res), msg=msg_)

    def assert_mock_args_in_last_call(self, mck, args=None, kwargs=None, msg=None):
        args_, kwargs_ = mck.call_args
        if args != None:
            for arg in args:
                msg_ = msg if msg != None else \
                'Argument %s was not in call to %s (called with args: %s)'%(arg, mck, args_)
                self.assertInList(arg, args_, msg=msg_)
        if kwargs != None:
            for key in list(kwargs.keys()):
                msg_ = msg if msg != None else \
                'Key Word Argument %s=%s was not in call to %s (called with args: %s)'% \
                (key, kwargs[key], mck, kwargs_)
                self.assertTrue(key in kwargs_, msg_)
                try:
                    self.assertEqual(kwargs[key], kwargs_[key], msg_)
                except Exception:
                    self.assertArrayAlmostEqual(kwargs[key], kwargs_[key], places=5)


class ModelTestCase(SEDTestCase):

    @classmethod
    def setUpClass(ModelTestCase):
        """ Setup the tmap grid as it is smaller and thus faster as kurucz"""
        model.set_defaults(grid='kurucztest')
        grid1 = dict(grid='tmaptest')
        grid2 = dict(grid='tmaptest')
        model.set_defaults_multiple(grid1,grid2)
        model.copy2scratch(z='*', Rv='*')

    def setUp(self):
        self.photbands = ['STROMGREN.U', '2MASS.H']

    def testGetItablePixSingle(self):
        """ model.get_itable_pix() single case """
        bgrid = {'teff': array([ 5500.,  6874.,  9645.,  7234., 5932.]),
                'logg': array([ 3.57,  4.00,  4.21,  3.84,  3.25]),
                'ebv': array([ 0.0018,  0.0077,  0.0112,  0.0046,  0.0110]),
                'rv': array([ 2.20, 2.40, 2.60, 2.80, 3.00]),
                'z': array([ -0.5, -0.4, -0.3, -0.2, 0.0])}

        flux_,Labs_ = model.get_itable_pix(photbands=self.photbands, **bgrid)

        flux = [[3255462., 13286738., 53641850., 16012786., 4652578.69189492],
                [967634., 1262321., 1789486., 1336016., 1066763.]]
        Labs = [0.819679, 1.997615, 7.739527, 2.450019, 1.107916]

        self.assertArrayAlmostEqual(flux_[0],flux[0],delta=100)
        self.assertArrayAlmostEqual(flux_[1],flux[1],delta=100)
        self.assertArrayAlmostEqual(Labs_,Labs,places=3)


    def testGetItablePixBinary(self):
        """ model.get_itable_pix() multiple case """
        bgrid = {'teff': array([ 22674.,  21774.,  22813.,  29343., 28170.]),
                'logg': array([ 5.75,  6.07,  6.03 ,  6.38,  5.97]),
                'ebv': array([ 0.0018,  0.0077,  0.0112,  0.0046,  0.0110]),
                'rad': array([ 4.96,  7.34,  4.56,  3.30,  8.54]),
                'teff2': array([ 38779.,  36447.,  32099. ,  31392., 35893.]),
                'logg2': array([ 4.67,  4.92,  5.46,  4.96,  4.85]),
                'ebv2': array([ 0.0018,  0.0077,  0.0112,  0.0046,  0.0110]),
                'rad2': array([ 6.99,  9.61,  9.55,  6.55,  3.99])}

        flux_,Labs_ = model.get_itable_pix(photbands=self.photbands, **bgrid)

        flux = [[1.69858e+11, 2.88402e+11, 1.98190e+11, 1.00088e+11, 1.39618e+11],
                [5.78019e+08, 1.04221e+09, 7.42561e+08, 3.63444e+08, 5.52846e+08]]
        Labs = [67046.5663, 115877.0362, 81801.3768, 39791.7467, 56369.3989]

        self.assertArrayAlmostEqual(flux_[0],flux[0],delta=0.01e+11)
        self.assertArrayAlmostEqual(flux_[1],flux[1],delta=0.01e+09)
        self.assertArrayAlmostEqual(Labs_,Labs,places=3)

    def testGetItableSingle(self):
        """ model.get_itable() single case """

        flux_,Labs_ = model.get_itable(photbands=self.photbands, teff=6874, logg=4.21,
                                                                   ebv=0.0077, z=-0.2)
        flux = [13428649.32576484, 1271090.21316342]
        Labs = 1.99865981674

        self.assertAlmostEqual(flux_[0],flux[0], delta=100)
        self.assertAlmostEqual(flux_[1],flux[1], delta=100)
        self.assertAlmostEqual(Labs_,Labs, delta=100)


    def testGetItableBinary(self):
        """ model.get_itable() multiple case """

        flux_,Labs_ = model.get_itable(photbands=self.photbands, teff=25000,
                        logg=5.12, ebv=0.001, teff2=33240, logg2=5.86, ebv2=0.001)

        flux = [3.31834622e+09, 1.25032866e+07]
        Labs = 1277.54498257

        self.assertAlmostEqual(flux_[0],flux[0], delta=100)
        self.assertAlmostEqual(flux_[1],flux[1], delta=100)
        self.assertAlmostEqual(Labs_,Labs, delta=100)

    def testGetTable(self):
        """ model.get_table() single case """

        wave, flux = model.get_table(teff=6874, logg=4.21, ebv=0.0)

        self.assertEqual(len(wave), len(flux))
        self.assertEqual(len(wave), 1221)
        self.assertAlmostEqual(wave[0], 90.9, delta=0.1)
        self.assertAlmostEqual(wave[-1], 1600000.0, delta=0.1)
        self.assertAlmostEqual(wave[400], 4370.0, delta=0.1)
        self.assertAlmostEqual(wave[800], 15925.0, delta=0.1)
        self.assertAlmostEqual(flux[0], 0.0, delta=0.0001)
        self.assertAlmostEqual(flux[-1], 0.020199, delta=0.0001)
        self.assertAlmostEqual(flux[400], 24828659.5845, delta=0.0001)
        self.assertAlmostEqual(flux[800], 1435461.60457, delta=0.0001)

    def testGetTableBinary(self):
        """ model.get_table() multiple case """

        wave,flux = model.get_table(teff=25000, logg=5.12, ebv=0.001,
                                    teff2=33240, logg2=5.86, ebv2=0.001)


        self.assertEqual(len(wave), len(flux))
        self.assertEqual(len(wave), 123104)
        self.assertAlmostEqual(wave[0], 1000, delta=0.1)
        self.assertAlmostEqual(wave[-1], 24999.8, delta=0.1)
        self.assertAlmostEqual(wave[40000], 8671.20117188, delta=0.001)
        self.assertAlmostEqual(wave[80000], 16509.4003906, delta=0.001)
        self.assertAlmostEqual(flux[0], 170415154318.0, delta=0.5)
        self.assertAlmostEqual(flux[-1], 2457826.26898, delta=0.0001)
        self.assertAlmostEqual(flux[40000], 141915936.111, delta=0.001)
        self.assertAlmostEqual(flux[80000], 12450102.801, delta=0.001)

class PixFitTestCase(SEDTestCase):

    @classmethod
    def setUpClass(PixFitTestCase):
        """ Setup the tmap grid as it is smaller and thus faster as kurucz"""
        model.set_defaults(grid='kurucztest')
        grid1 = dict(grid='kurucztest')
        grid2 = dict(grid='tmaptest')
        model.set_defaults_multiple(grid1,grid2)
        model.copy2scratch(z='*', Rv='*')


    def setUp(self):
        self.photbands = ['STROMGREN.U', '2MASS.H']

    def testGenerateGridPixSingle(self):
        """ fit.generate_grid_pix() single case """
        grid = fit.generate_grid_single_pix(self.photbands,teffrange=(5000,10000),
                    loggrange=(3.50,4.50),ebvrange=(0.0, 0.012),zrange=(-0.5,0.0),
                    rvrange=(-inf,inf),points=50)

        self.assertTrue('teff' in grid)
        self.assertTrue('logg' in grid)
        self.assertTrue('ebv' in grid)
        self.assertTrue('rv' in grid)
        self.assertTrue('z' in grid)

        self.assertFalse('teff2' in grid)
        self.assertFalse('rad' in grid)

        self.assertAllBetweenDiff(grid['teff'], 5000, 10000, places=-2)
        self.assertAllBetweenDiff(grid['logg'], 3.5, 4.5, places=1)
        self.assertAllBetweenDiff(grid['ebv'], 0.0, 0.012, places=3)
        self.assertAllBetweenDiff(grid['rv'], 2.1, 3.1, places=1)
        self.assertAllBetweenDiff(grid['z'], -0.5, 0.0, places=1)

    def testGenerateGridPixBinary(self):
        """ fit.generate_grid_pix() binary case """
        grid = fit.generate_grid_pix(self.photbands,teffrange=(5000,10000),
                    teff2range=(30000, 40000), loggrange=(3.5,4.50),
                    logg2range=(4.50, 6.50), ebvrange=(0.0, 0.012),
                    rvrange=(2.1,3.1), rv2range=(2.1,3.1),
                    masses=[0.47*1.63, 0.47], points=50)

        self.assertTrue('teff' in grid)
        self.assertTrue('logg' in grid)
        self.assertTrue('teff2' in grid)
        self.assertTrue('logg2' in grid)
        self.assertTrue('ebv' in grid)
        self.assertTrue('ebv2' in grid)
        self.assertTrue('rad' in grid)
        self.assertTrue('rad2' in grid)
        self.assertFalse('z' in grid)
        self.assertFalse('z2' in grid)

        self.assertListEqual(grid['ebv'].tolist(), grid['ebv2'].tolist(),
                       msg="ebv should be the same for both components.")
        self.assertListEqual(grid['rv'].tolist(), grid['rv2'].tolist(),
                      msg="Rv should be the same for both components.")

        self.assertAllBetweenDiff(grid['teff'], 5000, 10000, places=-2)
        self.assertAllBetweenDiff(grid['logg'], 3.5, 4.5, places=1)
        self.assertAllBetweenDiff(grid['ebv'], 0.0, 0.012, places=3)
        self.assertAllBetweenDiff(grid['rv'], 2.1, 3.1, places=1)
        self.assertAllBetweenDiff(grid['teff2'], 30000, 40000, places=-2)
        self.assertAllBetweenDiff(grid['logg2'], 4.5, 6.5, places=1)

        G, Msol, Rsol = constants.GG_cgs, constants.Msol_cgs, constants.Rsol_cgs
        rad = np.sqrt(G*0.47*1.63*Msol/10**grid['logg'])/Rsol
        rad2 = np.sqrt(G*0.47*Msol/10**grid['logg2'])/Rsol

        self.assertListEqual(rad.tolist(), grid['rad'].tolist())
        self.assertListEqual(rad2.tolist(), grid['rad2'].tolist())


    def testGenerateGridPixMultiple(self):
        """ fit.generate_grid_pix() multiple case """
        grid = fit.generate_grid_pix(self.photbands,teffrange=(5000,10000),
                    teff2range=(20000, 40000), loggrange=(3.5,4.5),radrange=(0.1,1.0),
                    logg2range=(4.50, 6.50), ebvrange=(0.0, 0.012),rad2range=(1.0,10.0),
                    points=50)

        self.assertTrue('teff' in grid)
        self.assertTrue('logg' in grid)
        self.assertTrue('teff2' in grid)
        self.assertTrue('logg2' in grid)
        self.assertTrue('ebv' in grid)
        self.assertTrue('ebv2' in grid)
        self.assertTrue('rad' in grid)
        self.assertTrue('rad2' in grid)
        self.assertFalse('z' in grid)
        self.assertFalse('z2' in grid)
        self.assertFalse('rv' in grid)
        self.assertFalse('rv2' in grid)

        self.assertListEqual(grid['ebv'].tolist(), grid['ebv2'].tolist())

        self.assertAllBetweenDiff(grid['teff'], 5000, 10000, places=-2)
        self.assertAllBetweenDiff(grid['logg'], 3.5, 4.5, places=1)
        self.assertAllBetweenDiff(grid['ebv'], 0.0, 0.012, places=3)
        self.assertAllBetweenDiff(grid['rad'], 0.1, 1.0, places=1)
        self.assertAllBetweenDiff(grid['teff2'], 20000, 40000, places=-2)
        self.assertAllBetweenDiff(grid['logg2'], 4.5, 6.5, places=1)
        self.assertAllBetweenDiff(grid['rad2'], 1.0, 10.0, places=0)

    @unittest.skipIf(noMock, "Mock not installed")
    def testiGridSearch(self):
        """ fit.igrid_search_pix() """
        meas = array([3.64007e-13, 2.49267e-13, 9.53516e-14] )
        emeas = array([3.64007e-14, 2.49267e-14, 9.53516e-15])
        photbands = ['STROMGREN.U', 'STROMGREN.B', 'STROMGREN.V']
        grid = {'teff': array([ 22674.,  21774.,  22813.,  29343., 28170.]),
                'logg': array([ 5.75,  6.07,  6.03 ,  6.38,  5.97]),
                'ebv': array([ 0.0018,  0.0077,  0.0112,  0.0046,  0.0110]),
                'z': array([0,0,0,0,0]),
                'rv': array([ 2.20, 2.40, 2.60, 2.80, 3.00])}
        syn_flux = array([[8.0218e+08, 7.2833e+08, 8.1801e+08, 1.6084e+09, 1.4178e+09],
                    [4.3229e+08, 4.0536e+08, 4.3823e+08, 7.0594e+08, 6.4405e+08],
                    [6.2270e+08, 5.7195e+08, 6.2482e+08, 1.0415e+09, 9.5594e+08]])
        lumis = array([232.5337, 200.0878, 238.7946, 625.3935, 533.8251])
        chisqs = array([260.1680, 255.4640, 251.1565, 221.6586, 233.4854])
        scales = array([3.9450e-22, 4.2714e-22, 3.8880e-22, 2.2365e-22, 2.4783e-22])
        e_scales = array([1.7789e-22, 1.9005e-22, 1.7449e-22, 1.0679e-22, 1.1745e-22])

        mock_model = self.create_patch(model, 'get_itable_pix', return_value=(syn_flux, lumis))
        mock_color = self.create_patch(filters, 'is_color', return_value=False)
        mock_stat = self.create_patch(fit, 'stat_chi2', return_value=(chisqs,scales,e_scales))

        chisqs_,scales_,e_scales_,lumis_ = fit.igrid_search_pix(meas, emeas, photbands, model_func=model.get_itable_pix, stat_func=fit.stat_chi2,  **grid)

        self.assertListEqual(chisqs.tolist(), chisqs_.tolist())
        self.assertListEqual(scales_.tolist(), scales.tolist())
        self.assertListEqual(e_scales_.tolist(), e_scales.tolist())
        self.assertListEqual(lumis_.tolist(), lumis.tolist())

        self.assert_mock_args_in_last_call(mock_model, kwargs=grid)
        self.assert_mock_args_in_last_call(mock_model, kwargs={'photbands':photbands})

        mock_stat.assert_called()

    def testCreateParameterDict(self):
        """ fit.create_parameter_dict() """

        pars = dict(teff_value=5000, teff_min=4000, teff_max=6000,
                    logg_value=4.5, logg_min=3.5, logg_max=5.5,
                    ebv_value=0.001, ebv_min=0.0, ebv_max=0.01, ebv_vary=False,
                    z_value=0, z_vary=False,
                    rad_expr='G*sqrt(2*logg)/m')

        parameters = fit.create_parameter_dict(**pars)

        self.assertTrue('value' in parameters)
        self.assertTrue('min' in parameters)
        self.assertTrue('max' in parameters)
        self.assertTrue('vary' in parameters)
        self.assertTrue('expr' in parameters)

        names = parameters['names']
        self.assertTrue('teff' in names)
        self.assertTrue('logg' in names)
        self.assertTrue('ebv' in names)
        self.assertTrue('z' in names)
        self.assertTrue('rad' in names)

        for key in list(parameters.keys()):
            self.assertTrue(type(parameters[key]) == np.ndarray)

        self.assertEqual(parameters['value'][names == 'teff'], 5000)
        self.assertEqual(parameters['value'][names == 'logg'], 4.5)
        self.assertEqual(parameters['value'][names == 'ebv'], 0.001)
        self.assertEqual(parameters['value'][names == 'z'], 0)
        self.assertEqual(parameters['min'][names == 'teff'], 4000)
        self.assertEqual(parameters['min'][names == 'logg'], 3.5)
        self.assertEqual(parameters['min'][names == 'ebv'], 0.0)
        self.assertEqual(parameters['max'][names == 'teff'], 6000)
        self.assertEqual(parameters['max'][names == 'logg'], 5.5)
        self.assertEqual(parameters['max'][names == 'ebv'], 0.01)
        self.assertEqual(parameters['vary'][names == 'ebv'], False)
        self.assertEqual(parameters['expr'][names == 'rad'], 'G*sqrt(2*logg)/m')




class MinimizeFitTestCase(SEDTestCase):

    def testCreateParameterDict(self):
        """ fit.create_parameter_dict() """
        pars = {'logg_value': 4.0, 'vrad_vary': False, 'z_value': -0.3, 'vrad_value': 0,
                'logg_vary': True, 'rv_min': 2.1, 'vrad_min': 0, 'vrad_max': 0, 'z_max': 0.0,
                'ebv_max': 0.015, 'teff_value': 6000, 'z_vary': True, 'rv_value': 2.4,
                'teff_vary': True, 'logg_min': 3.5, 'rv_max': 3.1, 'z_min': -0.5,
                'ebv_min': 0.005, 'teff_min': 5000, 'logg_max': 4.5, 'ebv_value': 0.007,
                'teff_max': 7000, 'rv_vary': True, 'ebv_vary': True}
        exp = {'max': array([3.1, 7000, 4.5, 0.015, 0, 0.0], dtype=object),
               'vary': array([True, True, True, True, False, True], dtype=object),
               'names': array(['rv', 'teff', 'logg', 'ebv', 'vrad', 'z'], dtype='|S4'),
               'value': array([2.4, 6000, 4.0, 0.007, 0, -0.3], dtype=object),
               'min': array([2.1, 5000, 3.5, 0.005, 0, -0.5], dtype=object)}

        res = fit.create_parameter_dict(**pars)

        self.assertArrayEqual(res['max'], exp['max'])
        self.assertArrayEqual(res['vary'], exp['vary'])
        self.assertArrayEqual(res['min'], exp['min'])
        self.assertArrayEqual(res['names'], exp['names'])
        self.assertArrayEqual(res['value'], exp['value'])

    @unittest.skipIf(noMock, "Mock not installed")
    def testiMinimize(self):
        """ fit.iminimize() normal mocked """

        gm_return = (['minimizer'], ['startpars'], ['newmodels'], ['chisqrs'])
        gifm_return = (['chisqr'], ['nfev'], ['scale'], ['labs'], ['grid'])
        cpd_return = dict(names=['teff'], value=[5000], min=[4000], max=[6000], vary=[True])

        mock_sfit_m = self.create_patch(sigproc.fit, 'minimize', return_value='minimizer')
        mock_sfit_gm = self.create_patch(sigproc.fit, 'grid_minimize', return_value=gm_return)
        mock_sfit_sp = self.create_patch(sigproc.fit.Function, 'setup_parameters')
        mock_fit_gifm = self.create_patch(fit, '_get_info_from_minimizer', return_value=gifm_return)
        mock_fit_cpd = self.create_patch(fit, 'create_parameter_dict', return_value=cpd_return)

        meas = array([1.25e9, 2.34e8])
        emeas = array([1.25e7, 2.34e6])
        photbands = array(['J', 'C'])

        #-- 1 point
        grid, chisqr, nfev, scale, lumis = fit.iminimize(meas,emeas,photbands, points=None, epsfcn=0.01 )

        self.assert_mock_args_in_last_call(mock_sfit_sp, kwargs=cpd_return)
        self.assert_mock_args_in_last_call(mock_sfit_m, args=[photbands, meas], kwargs=dict(epsfcn=0.01, weights=1/emeas))
        self.assert_mock_args_in_last_call(mock_fit_gifm, args=[['minimizer'], photbands, meas, emeas])
        self.assertListEqual(grid,['grid'])
        self.assertListEqual(chisqr,['chisqr'])
        self.assertListEqual(nfev,['nfev'])
        self.assertListEqual(scale,['scale'])
        self.assertListEqual(lumis,['labs'])

    @unittest.skipIf(noMock, "Mock not installed")
    def testiMinimizeGrid(self):
        """ fit.iminimize() grid mocked """
        gm_return = (['minimizer'], ['startpars'], ['newmodels'], ['chisqrs'])
        gifm_return = (['chisqr'], ['nfev'], ['scale'], ['labs'], ['grid'])
        cpd_return = dict(names=['teff'], value=[5000], min=[4000], max=[6000], vary=[True])

        mock_sfit_m = self.create_patch(sigproc.fit, 'minimize', return_value='minimizer')
        mock_sfit_gm = self.create_patch(sigproc.fit, 'grid_minimize', return_value=gm_return)
        mock_sfit_sp = self.create_patch(sigproc.fit.Function, 'setup_parameters')
        mock_fit_gifm = self.create_patch(fit, '_get_info_from_minimizer', return_value=gifm_return)
        mock_fit_cpd = self.create_patch(fit, 'create_parameter_dict', return_value=cpd_return)

        meas = array([1.25e9, 2.34e8])
        emeas = array([1.25e7, 2.34e6])
        photbands = array(['J', 'C'])

        #-- 10 point
        grid, chisqr, nfev, scale, lumis = fit.iminimize(meas,emeas,photbands, points=10, epsfcn=0.01 )

        self.assert_mock_args_in_last_call(mock_sfit_sp, kwargs=cpd_return)
        self.assert_mock_args_in_last_call(mock_sfit_gm, args=[photbands, meas], kwargs=dict(epsfcn=0.01, weights=1/emeas, points=10))
        self.assert_mock_args_in_last_call(mock_fit_gifm, args=[['minimizer'], photbands, meas, emeas])
        self.assertListEqual(grid,['grid'])
        self.assertListEqual(chisqr,['chisqr'])
        self.assertListEqual(nfev,['nfev'])
        self.assertListEqual(scale,['scale'])
        self.assertListEqual(lumis,['labs'])


class BuilderTestCase(SEDTestCase):

    @classmethod
    def setUpClass(BuilderTestCase):
        if not noMock:
            sesame.search = mock.Mock(return_value={'plx':(0.0,0.0)})
            builder.SED.load_photometry = mock.Mock(return_value=None)


    def setUp(self):
        self.sed = builder.SED(ID='TEST',load_fits=False)
        self.sed.master = {}
        self.sed.results = {'igrid_search':{}}

    def testInit(self):
        """ builder.sed.__init__() ToDo """
        self.assertTrue(self.sed.ID == 'TEST')

    def testCollectResults(self):
        """ builder.sed.collect_results() """
        bgrid = {'teff': array([ 22674.,  21774.,  22813.,  29343., 28170.]),
                 'logg': array([ 5.75, 6.07,  6.03,  6.38,  5.97]),
                 'ebv': array([ 0.0018, 0.0077,  0.0112,  0.0046,  0.0110]),
                 'rv': array([ 2.20, 2.40, 2.60, 2.80, 3.00]),
                 'z': array([0,0,0,0,0])}
        fgrid = {'chisq': array([1.,3.,2.,0.1,np.nan])}

        self.sed.collect_results(grid=bgrid, fitresults=fgrid, mtype='igrid_search', selfact='chisq')
        res = self.sed.results['igrid_search']['grid']

        self.assertListEqual(res['chisq'].tolist(), [3.0, 2.0, 1.0, 0.1])
        self.assertListEqual(res['teff'].tolist(), [21774.,22813.,22674.,29343.])
        self.assertListEqual(res['logg'].tolist(), [6.07,6.03,5.75,6.38])
        self.assertListEqual(res['ebv'].tolist(), [0.0077,0.0112,0.0018,0.0046])
        self.assertListEqual(res['rv'].tolist(), [2.40,2.60,2.20,2.80])
        self.assertListEqual(res['z'].tolist(), [0,0,0,0])
        self.assertListEqual(res['ci_raw'].tolist(), [0.,0.,0.,0.])
        self.assertListEqual(res['ci_red'].tolist(), [0.,0.,0.,0.])

    def testCalculateDegreesOfFreedom(self):
        """ builder.sed.calculateDF() """
        ranges = {'teffrange': (5000,6000), 'teff2range': (20000,50000), 'loggrange': (4.5,4.5),
                  'logg2range': (4.5,5.5), 'ebvrange': (0.0,0.01), 'ebv2range': (0.0,0.01),
                  'zrange': (-0.5,0.5)}

        df, dfinfo = self.sed.calculateDF(**ranges)

        self.assertEqual(df, 6)
        self.assertTrue('ebv' in dfinfo)
        self.assertFalse('ebv2' in dfinfo)

    @unittest.skipIf(noMock, "Mock not installed")
    def testCalculateStatistics(self):
        """ builder.sed.calculate_statistics()"""
        dtypes = [('teff','f8'), ('logg','f8'),('ebv','f8'),('rv','f8'),('z','f8'), ('chisq','f8')]
        grid = [array([ 22674.,  21774.,  22813.,  29343., 28170.]),
                array([ 5.75, 6.07,  6.03,  6.38,  5.97]),
                array([ 0.0018, 0.0077,  0.0112,  0.0046,  0.0110]),
                array([ 2.20, 2.40, 2.60, 2.80, 3.00]),
                array([0,0,0,0,0]),
                array([1.,3.,2.,0.1,10.0])]
        master = np.rec.fromarrays(grid,dtype=dtypes)
        master = mlab.rec_append_fields(master, 'ci_raw', np.zeros(len(master)))
        master = mlab.rec_append_fields(master, 'ci_red', np.zeros(len(master)))
        self.sed.results['igrid_search']['grid'] = master
        self.sed.master['include'] = [True,True,False,True,True]

        with mock.patch.object(builder.SED, 'calculateDF', return_value=5) as mock_method:
            self.sed.calculate_statistics(df=5)

        res = self.sed.results['igrid_search']
        raw = [0.6826894, 0.9167354, 0.8427007, 0.2481703, 0.9984345]
        red = [0.2481703, 0.4161175, 0.3452791, 0.0796556, 0.6826894]

        self.assertFalse(mock_method.called)
        self.assertArrayAlmostEqual(res['grid']['ci_raw'].tolist(), raw, places=5)
        self.assertArrayAlmostEqual(res['grid']['ci_red'].tolist(), red, places=5)
        self.assertEqual(res['factor'], 10.0)


    def testCalculateConfidenceIntervals(self):
        """ builder.sed.calculate_confidence_intervals()  ToDo """
        pass

    def testGenerateRanges(self):
        """ builder.sed.generate_ranges() ToDo """
        pass


    @unittest.skipIf(noMock, "Mock not installed")
    def testiGridSearch(self):
        """ builder.sed.igrid_search() mocked """

        ranges = {'teffrange':(20000,30000), 'loggrange':(5.5,6.5)}
        grid = {'teff': array([ 22674.,  21774.,  22813.,  29343., 28170.]),
                'logg': array([ 5.75, 6.07,  6.03,  6.38,  5.97])}
        fitres = [[1.,3.,2.,0.1,10.0],[1.,3.,2.,0.1,10.0],[1.,3.,2.,0.1,10.0],[1.,3.,2.,0.1,10.0]]
        ci = dict(name=['teff', 'logg'], value=[30000, 5.5], cilow=[29000, 5.0], cihigh=[31000, 6.0])

        mock_sed_gr = self.create_patch(builder.SED, 'generate_ranges', return_value=ranges)
        mock_fit_ggp = self.create_patch(fit, 'generate_grid_pix', return_value=grid)
        mock_fit_isp = self.create_patch(fit, 'igrid_search_pix', return_value=fitres)
        mock_sed_cr = self.create_patch(builder.SED, 'collect_results')
        mock_sed_cs = self.create_patch(builder.SED, 'calculate_statistics')
        mock_sed_cci = self.create_patch(builder.SED, 'calculate_confidence_intervals', return_value=ci)
        mock_sed_sci = self.create_patch(builder.SED, 'store_confidence_intervals')
        mock_sed_sbm = self.create_patch(builder.SED, 'set_best_model')
        mock_p2s = self.create_patch(builder, 'photometry2str', return_value='TEST')

        self.sed.master = {'include':0, 0:0, 'cmeas':[0,0], 'e_cmeas':[0,0], 'photband':[0,0]}
        self.sed.igrid_search(teffrange=(10,20), loggrange=(4.5,4.5), ebvrange=(0,0.1), zrange=(0,0),rvrange=(3.1,3.1),vradrange=(0,0), df=4, set_model=True)

        mock_sed_gr.assert_called_with(teffrange=(10,20), loggrange=(4.5,4.5), ebvrange=(0,0.1), zrange=(0,0),rvrange=(3.1,3.1),vradrange=(0,0))
        mock_fit_ggp.assert_called()
        mock_fit_isp.assert_called()
        mock_sed_cr.assert_called()
        self.assert_mock_args_in_last_call(mock_sed_cr, kwargs={'grid':grid})
        mock_sed_cs.assert_called()
        self.assert_mock_args_in_last_call(mock_sed_cs, kwargs={'df':4})
        mock_sed_sci.assert_called()
        self.assert_mock_args_in_last_call(mock_sed_sci, kwargs=ci)
        mock_sed_cci.assert_called()
        mock_sed_sbm.assert_called()


class XIntegrationTestCase(SEDTestCase):

    photbands = ['STROMGREN.U', 'STROMGREN.B', 'STROMGREN.V', 'STROMGREN.Y',
                 '2MASS.H', '2MASS.J', '2MASS.KS']
    measHot = None
    measCold = None
    measBin = None

    @classmethod
    def setUpClass(cls):
        if not noMock:
            sesame.search = mock.Mock(return_value={'plx':(0.0,0.0)})
        if not noIntegration:
            # ==== COLD model ====
            model.set_defaults(grid='kurucztest')
            model.copy2scratch(z='*', Rv='*')
            measCold = model.get_itable_pix(photbands=cls.photbands, teff=array([6000]), \
                                logg=array([4.0]),ebv=array([0.01]), rv=array([2.8]), z=array([-0.25]))[0][:,0]

            np.random.seed(111)
            cls.measCold = measCold

            # ==== HOT model ====
            model.set_defaults(grid='tmaptest')
            model.copy2scratch(z='*', Rv='*')
            measHot = model.get_itable_pix(photbands=cls.photbands, teff=array([30000]), \
                                logg=array([5.5]),ebv=array([0.01]), rv=3.1, z=0.0)[0][:,0]

            np.random.seed(111)
            cls.measHot = measHot

            # ==== BINARY model ====
            grid1 = dict(grid='kurucztest')
            grid2 = dict(grid='tmaptest')
            model.set_defaults_multiple(grid1,grid2)
            model.clean_scratch(z='*', Rv='*')
            model.copy2scratch(z='*', Rv='*')

            G, Msol, Rsol = constants.GG_cgs, constants.Msol_cgs, constants.Rsol_cgs
            masses = [0.85, 0.50]
            rad = array([np.sqrt(G*masses[0]*Msol/10**4.0)/Rsol])
            rad2 = array([np.sqrt(G*masses[1]*Msol/10**5.5)/Rsol])

            measBin = model.get_itable_pix(photbands=cls.photbands, teff=array([6000]), \
                                logg=array([4.0]),ebv=array([0.01]), teff2=array([30000]),
                                logg2=array([5.5]), ebv2=array([0.01]), rad=rad,
                                rad2=rad2)[0][:,0]

            np.random.seed(111)
            cls.measBin = measBin
            cls.masses = masses

    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiGrid_searchSingleCold(self):
        """ INTEGRATION igrid_search single star (kurucz)"""
        sed = builder.SED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measCold + np.random.uniform(0, 0.01, size=len(self.measCold)) * self.measCold
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)

        model.set_defaults(grid='kurucztest')
        model.copy2scratch(z='*', Rv='*')

        np.random.seed(111)
        sed.igrid_search(points=100000,teffrange=(5000, 7000),loggrange=(3.5, 4.5),
                         ebvrange=(0.005, 0.015),zrange=(-0.5,0.0),rvrange=(2.1,3.1),
                         vradrange=(0,0),df=None,CI_limit=0.95,set_model=True)

        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff'], 6000, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg'], 3.98, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv'], 0.011, delta=0.02)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['rv'], 2.13, delta=0.5)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['z'], -0.28, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_l'], 5949, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_l'], 3.62, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_l'], 0.005, delta=0.02)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['rv_l'], 2.1, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['z_l'], -0.46, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_u'], 6060, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_u'], 4.5, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_u'], 0.015, delta=0.02)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['rv_u'], 3.1, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['z_u'], -0.05, delta=0.1)

        # check that the best model is stored
        self.assertTrue('model' in sed.results['igrid_search'])
        self.assertTrue('synflux' in sed.results['igrid_search'])
        self.assertTrue('chi2' in sed.results['igrid_search'])
        self.assertEqual(len(sed.results['igrid_search']['model']), 3, msg='stored model has wrong number of collumns (should be 3)')
        self.assertEqual(len(sed.results['igrid_search']['synflux']), 3, msg='stored synflux has wrong number of collumns (should be 3)')


    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiGrid_searchSingleHot(self):
        """ INTEGRATION igrid_search single star (tmap) """
        sed = builder.SED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measHot + np.random.uniform(0, 0.04, size=len(self.measHot)) * self.measHot
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)

        model.set_defaults(grid='tmaptest')
        model.copy2scratch(z='*', Rv='*')

        np.random.seed(111)
        sed.igrid_search(points=100000,teffrange=(25000, 35000),loggrange=(5.0, 6.0),
                         ebvrange=(0.005, 0.015),zrange=(0,0),rvrange=(3.1,3.1),
                         vradrange=(0,0),df=None,CI_limit=0.95,set_model=True)

        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff'], 30200, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg'], 5.67, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv'], 0.0078, delta=0.02)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_l'], 29337, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_l'], 5.0, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_l'], 0.005, delta=0.02)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_u'], 31623, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_u'], 6.0, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_u'], 0.015, delta=0.02)

        # check that the best model is stored
        self.assertTrue('model' in sed.results['igrid_search'])
        self.assertTrue('synflux' in sed.results['igrid_search'])
        self.assertTrue('chi2' in sed.results['igrid_search'])
        self.assertEqual(len(sed.results['igrid_search']['model']), 3, msg='stored model has wrong number of collumns (should be 3)')
        self.assertEqual(len(sed.results['igrid_search']['synflux']), 3, msg='stored synflux has wrong number of collumns (should be 3)')

    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiGrid_searchBinary(self):
        """ INTEGRATION igrid_search binary star (kurucz-tmap) """
        sed = builder.BinarySED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measBin + np.random.uniform(0, 0.04, size=len(self.measBin)) * self.measBin
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)

        grid1 = dict(grid='kurucztest')
        grid2 = dict(grid='tmaptest')
        model.set_defaults_multiple(grid1,grid2)

        np.random.seed(111)
        sed.igrid_search(points=100000,teffrange=[(5000, 7000),(25000, 35000)],
                         loggrange=[(3.5, 4.5),(5.0, 6.0)],
                         ebvrange=[(0.005, 0.015), (0.005, 0.015)],
                         zrange=[(0,0),(0,0)],
                         radrange=[(0,10),(0,10)],
                         masses=self.masses,
                         df=None,CI_limit=0.95,set_model=False)

        teff = sed.results['igrid_search']['CI']['teff']
        logg = sed.results['igrid_search']['CI']['logg']
        ebv = sed.results['igrid_search']['CI']['ebv']
        teff2 = sed.results['igrid_search']['CI']['teff2']
        logg2 = sed.results['igrid_search']['CI']['logg2']
        ebv2 = sed.results['igrid_search']['CI']['ebv2']

        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff'], 6134, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg'], 4.38, delta=0.25)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv'], 0.008, delta=0.002)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_l'], 5826, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_l'], 4.35, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_l'], 0.005, delta=0.002)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff_u'], 6291, delta=50)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg_u'], 4.40, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv_u'], 0.015, delta=0.002)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff2'], 32195, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg2'], 5.97, delta=0.25)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv2'], 0.008, delta=0.002)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff2_l'], 26820, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg2_l'], 5.70, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv2_l'], 0.005, delta=0.002)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['teff2_u'], 33025, delta=250)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['logg2_u'], 5.99, delta=0.1)
        self.assertAlmostEqual(sed.results['igrid_search']['CI']['ebv2_u'], 0.015, delta=0.002)

    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiMinimizeSingleCold(self):
        """ INTEGRATION iminimize single star (kurucz) """
        sed = builder.SED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measCold + np.random.uniform(0, 0.04, size=len(self.measCold)) * self.measCold
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)

        model.set_defaults(grid='kurucztest')
        model.copy2scratch(z='*', Rv='*')

        np.random.seed(111)
        sed.iminimize(teff=6000, logg=4.0, ebv=0.007, z=-0.3, rv=2.4, vrad=0,
                      teffrange=(5000, 7000),loggrange=(3.5, 4.5),zrange=(-0.5,0.0),
                      ebvrange=(0.005, 0.015), rvrange=(2.1,3.1),vradrange=(0,0),
                      points=None,df=None,CI_limit=0.60,calc_ci=True, set_model=True)

        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff'], 6036, delta=50)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['logg'], 4.19, delta=0.1)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['ebv'], 0.015, delta=0.02)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['z'], -0.21, delta=0.1)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['rv'], 2.1, delta=0.3)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['scale'], 1, delta=0.5)

        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff_l'], 6025, delta=50)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff_u'], 6036, delta=50)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['scale_l'], 1, delta=0.5)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['scale_u'], 1, delta=0.5)

        self.assertEqual(sed.results['iminimize']['grid']['teffstart'][0], 6000)
        self.assertEqual(sed.results['iminimize']['grid']['loggstart'][0], 4.0)
        self.assertEqual(sed.results['iminimize']['grid']['ebvstart'][0], 0.007)
        self.assertEqual(sed.results['iminimize']['grid']['zstart'][0], -0.3)
        self.assertEqual(sed.results['iminimize']['grid']['rvstart'][0], 2.4)
        self.assertAlmostEqual(sed.results['iminimize']['grid']['chisq'][0], 3.9, delta=1)

        self.assertTrue('model' in sed.results['iminimize'])
        self.assertTrue('synflux' in sed.results['iminimize'])
        self.assertTrue('chi2' in sed.results['iminimize'])
        self.assertEqual(len(sed.results['iminimize']['model']), 3, msg='stored model has wrong number of collumns (should be 3)')
        self.assertEqual(len(sed.results['iminimize']['synflux']), 3, msg='stored synflux has wrong number of collumns (should be 3)')


    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiMinimizeSingleHot(self):
        """ INTEGRATION iminimize single star (tmap) """
        sed = builder.SED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measHot + np.random.uniform(0, 0.04, size=len(self.measHot)) * self.measHot
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)

        model.set_defaults(grid='tmaptest')
        model.copy2scratch(z='*', Rv='*')

        np.random.seed(111)
        sed.iminimize(teff=27000, logg=5.1, ebv=0.01, z=0, rv=3.1, vrad=0,
                      teffrange=(25000, 35000),loggrange=(5.0, 6.0),
                      ebvrange=(0.005, 0.015),zrange=(0,0),rvrange=(3.1,3.1),
                      vradrange=(0,0),df=None,CI_limit=0.95,set_model=False)

        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff'], 30250, delta=100)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['logg'], 5.66, delta=0.1)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['ebv'], 0.008, delta=0.02)
        self.assertAlmostEqual(sed.results['iminimize']['grid']['chisq'][0], 3.8, delta=1)


    @unittest.skipIf(noIntegration, "Integration tests are skipped.")
    def testiMinimizeBinary(self):
        """ INTEGRATION iminimize binary star (kurucz-tmap) """
        sed = builder.BinarySED(ID='TEST', load_fits=False)
        np.random.seed(111)
        meas = self.measBin + np.random.uniform(0, 0.01, size=len(self.measBin)) * self.measBin
        emeas = meas / 100.0
        units = ['erg/s/cm2/AA' for i in meas]
        source = ['SYNTH' for i in meas]
        sed.add_photometry_fromarrays(meas, emeas, units, self.photbands, source)
        sed.set_photometry_scheme('abs')

        np.random.seed(111)
        sed.iminimize(teff=(5300,29000), logg=(4.0,5.3), ebv=(0.01,0.01), z=(0,0), rv=(3.1,3.1),
                    vrad=(0,0),teffrange=[(5000,7000),(25000, 35000)],loggrange=[(3.5,4.5),
                    (5.0, 6.0)], ebvrange=[(0.01,0.01),(0.01, 0.01)] ,zrange=[(0,0),(0,0)],
                    rvrange=[(3.1,3.1),(3.1,3.1)], vradrange=[(0,0),(0,0)], df=2,
                    CI_limit=None, masses = self.masses,calc_ci=False, set_model=True)

        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff'], 6023, delta=100)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['logg'], 3.95, delta=0.1)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['ebv'], 0.01, delta=0.02)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['teff2'], 30506, delta=100)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['logg2'], 5.47, delta=0.1)
        self.assertAlmostEqual(sed.results['iminimize']['CI']['scale'], 0.9, delta=0.5)

        self.assertTrue('model' in sed.results['iminimize'])
        self.assertTrue('synflux' in sed.results['iminimize'])
        self.assertTrue('chi2' in sed.results['iminimize'])
        self.assertEqual(len(sed.results['iminimize']['model']), 3, msg='stored model has wrong number of collumns (should be 3)')
        self.assertEqual(len(sed.results['iminimize']['synflux']), 3, msg='stored synflux has wrong number of collumns (should be 3)')



















