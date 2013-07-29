import numpy as np
import pylab as pl
from ivs.spectra import tools

import unittest

class SpectrumTestCase(unittest.TestCase):
    """Add some extra usefull assertion methods to the testcase class"""
    
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

class MergeCosmicClippingTestCase(SpectrumTestCase):
    
    def setUp(self):
        np.random.seed(1111)
        num = 1000
        level = np.array([134, 53, 261, 111])
        noise = 0.01 * level
        
        # create standard noise
        wave = np.linspace(100,200, num=num)
        flux1 = np.random.normal(level[0], noise[0], num)
        flux2 = np.random.normal(level[1], noise[1], num)
        flux3 = np.random.normal(level[2], noise[2], num)
        flux4 = np.random.normal(level[3], noise[3], num)
        
        # create cosmics
        flux1[467] *= 2.5
        flux2[197] *= 2.6
        flux3[816] *= 2.4
        flux4[23] *= 2.8
        
        # create spectral feature
        flux1[365] *= 1.40
        flux2[365] *= 1.43
        flux3[365] *= 1.38
        flux4[365] *= 1.39
        
        self.waves = np.array([wave, wave, wave, wave])
        self.fluxes = np.array([flux1, flux2, flux3, flux4])
        self.num = num
        self.avg = np.sum(level) * 1.05
        
    
    def testMergeSpectra(self):
        """ spectra.tools.merge_cosmic_clipping() no vrad """
        
        wave, flux = tools.merge_cosmic_clipping(self.waves, self.fluxes, sigma=2.0)
        avg = self.avg
        
        #-- Test removed cosmics
        self.assertTrue(flux[467] < avg, 'Cosmic not removed')
        self.assertTrue(flux[197] < avg, 'Cosmic not removed')
        self.assertTrue(flux[816] < avg, 'Cosmic not removed')
        self.assertTrue(flux[23] < avg, 'Cosmic not removed')
        
        #-- Test kept spectral feature
        self.assertTrue(flux[365] > 1.3 * avg and flux[365] < 1.5 * avg, 'Feature wrongly removed')
        
    def testMergeSpectraVrad(self):
        """ spectra.tools.merge_cosmic_clipping() with vrad """
        
        pl.plot(self.waves[0], np.sum(self.fluxes, axis=0))
        vrads = [555., 565., 545., 560.]
        wave, flux = tools.merge_cosmic_clipping(self.waves, self.fluxes, vrads=vrads, sigma=2.0)
        avg = self.avg
        
        #-- Test removed cosmics
        self.assertTrue(flux[467] < avg, 'Cosmic not removed')
        self.assertTrue(flux[197] < avg, 'Cosmic not removed')
        self.assertTrue(flux[816] < avg, 'Cosmic not removed')
        self.assertTrue(flux[23] < avg, 'Cosmic not removed')
        
        #-- Test kept spectral feature
        self.assertTrue(flux[365] > 1.3 * avg and flux[365] < 1.5 * avg, 'Feature wrongly removed')
        
        #-- Test shifted spectra
        self.assertAlmostEqual(wave[0], 99.8148719272, places=3, msg='Spectrum not shifted')
        self.assertAlmostEqual(wave[-1], 199.629743854, places=3, msg='Spectrum not shifted')
        f = flux[(wave > 136.2) & (wave < 136.36)]
        self.assertTrue(f[0] > 1.3 * avg, 'Spectrum not shifted')

    def testMergeSpectraFO(self):
        """ spectra.tools.merge_cosmic_clipping() full output """
        
        wave, flux, accepted, rejected = tools.merge_cosmic_clipping(self.waves, self.fluxes,
                                                                 sigma=2.0, full_output=True)
        
        #-- Test Rejected
        rej = np.array([467, 197, 816,  23])
        self.assertArrayEqual(rejected[0], np.arange(len(self.waves)), msg='Rejected not oke')
        self.assertArrayEqual(rejected[1], rej, msg='Rejected not oke')
        
        #-- Test accepted
        #print accepted
        self.assertEqual(len(accepted[1][accepted[0] == 0]), 999, msg='Accepted 0 not oke')
        self.assertEqual(len(accepted[1][accepted[0] == 1]), 999, msg='Accepted 1 not oke')
        self.assertEqual(len(accepted[1][accepted[0] == 2]), 999, msg='Accepted 2 not oke')
        self.assertEqual(len(accepted[1][accepted[0] == 3]), 999, msg='Accepted 3 not oke')
        a = list(np.arange(self.num))
        a.remove(467)
        self.assertArrayEqual(accepted[1][accepted[0] == 0], a, msg='Accepted 0 not oke')
        a = list(np.arange(self.num))
        a.remove(197)
        self.assertArrayEqual(accepted[1][accepted[0] == 1], a, msg='Accepted 1 not oke')
        a = list(np.arange(self.num))
        a.remove(816)
        self.assertArrayEqual(accepted[1][accepted[0] == 2], a, msg='Accepted 2 not oke')
        a = list(np.arange(self.num))
        a.remove(23)
        self.assertArrayEqual(accepted[1][accepted[0] == 3], a, msg='Accepted 3 not oke')
        
        