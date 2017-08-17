import numpy as np
import pylab as pl
from ivs.inout import fits
from ivs.spectra import tools

import unittest

class SpectrumTestCase(unittest.TestCase):
    """Add some extra usefull assertion methods to the testcase class"""
    
    def assertNoCosmic(self, wave, flux, wrange, sigma=3, msg=None):
        flux = flux[(wave>=wrange[0]-0.5) & (wave<=wrange[1]+0.5)]
        wave = wave[(wave>=wrange[0]-0.5) & (wave<=wrange[1]+0.5)]
        avg = np.median(flux)
        std = np.std(flux)
        
        msg_ = "Cosmic detected in %s <-> %s"%(wrange[0], wrange[1])
        if msg != None: msg_ = msg_ + ", " + str(msg)
        self.assertTrue(np.all(flux < avg + sigma * std), msg=msg_)
    
    def assertAbsorptionLine(self, wave, flux, line, cont, depth=0.75, msg=None):
        lf = np.min(flux[(wave>=line[0]) & (wave<=line[1])])
        cf = np.average(flux[(wave>=cont[0]) & (wave<=cont[1])])
        
        msg_ = "Absorption line missing in %s <-> %s"%(line[0], line[1])
        if msg != None: msg_ = msg_ + ", " + str(msg)
        self.assertTrue( lf <= cf * depth, msg=msg_ )
    
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

class MergeCosmicClippingTestCase1(SpectrumTestCase):
    """ Testcase using a synthetic spectrum with 1 emission line and 4 cosmics """
    
    @classmethod  
    def setUpClass(cls):
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
        
        waves = np.array([wave, wave, wave, wave])
        fluxes = np.array([flux1, flux2, flux3, flux4])
        vrads = [555., 565., 545., 560.]
        
        wave1, flux1, accepted, rejected = tools.merge_cosmic_clipping(waves, fluxes,
                                     sigma=2.0, base='median', offset='std', window=51,
                                     runs=2, full_output=True)
        wave2, flux2 = tools.merge_cosmic_clipping(waves, fluxes, vrads=vrads, 
                            sigma=2.0, base='median', offset='std', window=51, runs=2)
        
        cls.wave1 = wave1
        cls.flux1 = flux1
        cls.wave2 = wave2
        cls.flux2 = flux2
        cls.accepted = accepted
        cls.rejected = rejected
        cls.waves = waves
        cls.fluxes = fluxes
        cls.vrads = vrads
        cls.num = num
        cls.avg = np.sum(level) * 1.05
    
    def testRemoveCosmics(self):
        """ spectra.tools.merge_cosmic_clipping() synth remove cosmics """
        self.assertTrue(self.flux1[467] < self.avg, 'Cosmic NOT removed at 467')
        self.assertTrue(self.flux1[197] < self.avg, 'Cosmic NOT removed at 197')
        self.assertTrue(self.flux1[816] < self.avg, 'Cosmic NOT removed at 816')
        self.assertTrue(self.flux1[23] < self.avg, 'Cosmic NOT removed at 23')
    
    def testKeepSpectralFeatures(self):
        """ spectra.tools.merge_cosmic_clipping() synth keep features """
        self.assertTrue(self.flux1[365] > 1.3 * self.avg and self.flux1[365] < 1.5 * self.avg,
                        'Feature wrongly removed')
    
    def testRemoveCosmicsVrad(self):
        """ spectra.tools.merge_cosmic_clipping() synth Vrad remove cosmics """
        self.assertTrue(self.flux2[467] < self.avg, 'Cosmic NOT removed at 467')
        self.assertTrue(self.flux2[197] < self.avg, 'Cosmic NOT removed at 197')
        self.assertTrue(self.flux2[816] < self.avg, 'Cosmic NOT removed at 816')
        self.assertTrue(self.flux2[23] < self.avg, 'Cosmic NOT removed at 23')
    
    def testKeepSpectralFeaturesVrad(self):
        """ spectra.tools.merge_cosmic_clipping() synth Vrad keep features """
        self.assertTrue(self.flux2[365] > 1.3 * self.avg and self.flux2[365] < 1.5 * self.avg, 
                        'Feature wrongly removed')
    
    def testUseWaveFirstSpectrum(self):
        """ spectra.tools.merge_cosmic_clipping() synth return wavescale first spectrum """
        self.assertArrayAlmostEqual(self.wave1, self.waves[0], places=3, 
                                    msg='Did not return wavelength of first spectrum')
        
    def testShiftSpectraVrad(self):
        """ spectra.tools.merge_cosmic_clipping() synth shift spectra with Vrad """
        self.assertAlmostEqual(self.wave2[0], 99.8148719272, places=3, msg='Spectrum not shifted')
        self.assertAlmostEqual(self.wave2[-1], 199.629743854, places=3, msg='Spectrum not shifted')
        f = self.flux2[(self.wave2 > 136.2) & (self.wave2 < 136.36)]
        self.assertTrue(f[0] > 1.3 * self.avg, 'Spectrum not shifted')
    
    def testFullOutput(self):
        """ spectra.tools.merge_cosmic_clipping() synth Full Output """
        rejected = self.rejected
        accepted = self.accepted
        #-- Test Rejected
        rej = np.array([467, 197, 816,  23])
        self.assertArrayEqual(rejected[0], np.arange(len(self.waves)), msg='Rejected not oke')
        self.assertArrayEqual(rejected[1], rej, msg='Rejected not oke')
        
        #-- Test accepted
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

class MergeCosmicClippingTestCase2(SpectrumTestCase):
    """ Testcase using real spectra of Feige 66 (5000 - 7000 AA) """
    
    @classmethod  
    def setUpClass(cls):
        temp = '/STER/mercator/hermes/%s/reduced/%s_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'
        objlist = [('20090619','237033'), ('20090701','240226'),
                    ('20090712','241334'), ('20090712','241335'),
                    ('20100107','00268012'), ('20100120','00272619'),
                    ('20100120','00272620'), ('20100203','00273577'),
                    ('20100203','00273578'), ('20100303','00275671'),
                    ('20100303','00275672'), ('20100410','00281505'),
                    ('20100519','00284636'), ('20110222','00334558'),
                    ('20110319','00336547'), ('20110324','00339848'),
                    ('20110401','00342273'), ('20110402','00342363'),
                    ('20110406','00342699'), ('20110408','00342896'),
                    ('20120107','00391289'), ('20120110','00391633'),
                    ('20120116','00392217'), ('20120127','00393151'),
                    ('20120209','00394175'), ('20120330','00399697'),
                    ('20120420','00404769'), ('20120506','00406531'),
                    ('20130106','00445346'), ('20130215','00452556'),
                    ('20130406','00457718'), ('20130530','00474128')]
        mergeList = [temp%o for o in objlist]
        
        waves, fluxes = [], []
        for ifile in mergeList:
            w, f = fits.read_spectrum(ifile)
            f = f[(w>5000) & (w<7000)]
            w = w[(w>5000) & (w<7000)]
            waves.append(w)
            fluxes.append(f)
         
        wave, flux, accepted, rejected = tools.merge_cosmic_clipping(waves, fluxes, sigma=5.0,
                            base='average', offset='std', window=51, runs=2, full_output=True)
        
        cls.nspek = len(objlist)
        cls.wave = wave
        cls.flux = flux
        cls.accepted = accepted
        cls.rejected = rejected
        
    def testRemoveCosmics(self):
        """ spectra.tools.merge_cosmic_clipping() Feige 66 remove cosmics """
        self.assertNoCosmic(self.wave, self.flux, (6161.7, 6161.9))
        self.assertNoCosmic(self.wave, self.flux, (6166.5, 6166.8), sigma=4)
        self.assertNoCosmic(self.wave, self.flux, (6200.2, 6200.5))
        self.assertNoCosmic(self.wave, self.flux, (6546.9, 6547.1))
    
    def testKeepSpectralFeatures(self):
        """ spectra.tools.merge_cosmic_clipping() Feige 66 keep features """
        self.assertAbsorptionLine(self.wave, self.flux, (5159.8, 5160.3), (5157.8, 5159.4),
                                  depth=0.8)
        self.assertAbsorptionLine(self.wave, self.flux, (5874.9, 5876.2), (5870.7, 5872.9),
                                  depth=0.6)
        self.assertAbsorptionLine(self.wave, self.flux, (6001.7, 6002.0), (6000., 6001.45),
                                  depth=0.95)
    def testFullOutput(self):
        """ spectra.tools.merge_cosmic_clipping() Feige 66 Full Output """
        rejected = self.rejected
        accepted = self.accepted
        
        self.assertTrue(len(accepted[0]) > 2070043-100 and len(accepted[0]) < 2070043+100, 
                        msg=' Exceptional deviation from expected acceptions! ')
        self.assertTrue(len(rejected[0]) > 250 and len(rejected[0]) < 400, 
                        msg=' Exceptional deviation from expected rejections! ')
