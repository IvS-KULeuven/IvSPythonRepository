import os
import nose
import shutil
import getpass
import numpy as np
import pylab as pl
from ivs.catalogs import hermes
from ivs.io import fits

import unittest
try:
    import mock
    from mock import patch
    noMock = False
except Exception:
    noMock = True


class HermesTestCase(unittest.TestCase):
    """Add some extra usefull assertion methods to the testcase class"""
    
    test_data_dir = '/home/jorisv/IVS_python/ivs_test_data/'
    testdirectories = []
    
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
    
    @nose.tools.nottest
    def make_testDir(self):
        
        testdir = '/home/{:}/{:.0f}/'.format(getpass.getuser(), np.random.uniform()*10**6)
        while os.path.exists(testdir):
            testdir = '/home/{:}/{:.0f}/'.format(getpass.getuser(), np.random.uniform()*10**6)
            
        os.makedirs(testdir)
        
        self.testdirectories.append(testdir)
        
        return testdir
    
    @nose.tools.nottest
    def delete_testDir(self):
        
        for dirname in self.testdirectories:
            os.rmdir(dirname)
                
class TestCase1MergeHermesSpectra(HermesTestCase):
    
    def testOrderMerged(self):
        """Merge order-merged hermes spectra NO vrad (Feige 66)"""
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
        
        #wtotal, f = fits.read_spectrum(mergeList[0])
        #ftotal = np.zeros_like(f)
        #for ifile in mergeList:
            #w, f = fits.read_spectrum(ifile)
            #ftotal += np.interp(wtotal,w,f) 
        #pl.plot(wtotal, ftotal)
        
        wave, flux, header = hermes.merge_hermes_spectra(mergeList, sigma=3, offset=None,
                                                         base='average')
        #pl.plot(wave, flux)
        #pl.show()
        
        #-- Test that cosmics are removed
        self.assertNoCosmic(wave, flux, (4362.7, 4362.9))
        self.assertNoCosmic(wave, flux, (4974.2, 4974.4))
        self.assertNoCosmic(wave, flux, (5656.1, 5656.3))
        self.assertNoCosmic(wave, flux, (6161.7, 6161.9))
        self.assertNoCosmic(wave, flux, (6166.5, 6166.8), sigma=4)
        self.assertNoCosmic(wave, flux, (6200.2, 6200.5))
        self.assertNoCosmic(wave, flux, (6546.9, 6547.1))
        self.assertNoCosmic(wave, flux, (7996.8, 7997.0))
        self.assertNoCosmic(wave, flux, (8001.7, 8001.9))
        
        
        #-- Test that absorption lines are still there
        self.assertAbsorptionLine(wave, flux, (3933.4, 3933.7), (3934.1, 3934.9), depth=0.6)
        self.assertAbsorptionLine(wave, flux, (5159.8, 5160.3), (5157.8, 5159.4), depth=0.8)
        self.assertAbsorptionLine(wave, flux, (5874.9, 5876.2), (5870.7, 5872.9), depth=0.6)
        self.assertAbsorptionLine(wave, flux, (6001.7, 6002.0), (6000., 6001.45), depth=0.95)
        
    
    def testOrderMergedVrad(self):
        """Merge order-merged hermes spectra with vrad (KIC 9540226)"""
        
        pre = '/STER/mercator/hermes/'
        post_merged = '_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'
        #
        
        objlist = ['20120730/reduced/00415609',
                    '20120808/reduced/00416266',
                    '20120808/reduced/00416267',
                    '20120818/reduced/00417309',
                    '20120818/reduced/00417310',
                    '20120905/reduced/00419141',
                    '20120911/reduced/00419682',
                    '20120913/reduced/00419913',
                    '20120913/reduced/00419914',
                    '20120924/reduced/00421148']
        mergeList = [pre+o+post_merged for o in objlist]
        
        vrad = [-21.190, -17.473, -17.528, -11.621, -11.613, 5.974, 13.297, 15.793, 15.810, 18.739]
        
        #wtotal, f = fits.read_spectrum(mergeList[0])
        #ftotal = np.zeros_like(f)
        #for ifile in mergeList:
            #w, f = fits.read_spectrum(ifile)
            #ftotal += np.interp(wtotal,w,f) 
        #pl.plot(wtotal, ftotal)
        
        wave, flux, header = hermes.merge_hermes_spectra(mergeList, vrads=vrad, sigma=1.5,
                                                         offset='none', base='average', runs=3)
        #pl.plot(wave, flux)
        #pl.show()
        
        #-- Test that cosmics are removed
        self.assertNoCosmic(wave, flux, (4060.1, 4060.3))
        self.assertNoCosmic(wave, flux, (4212.3, 4212.4))
        self.assertNoCosmic(wave, flux, (4320.0, 4320.2), sigma=2)
        self.assertNoCosmic(wave, flux, (7056.5, 7056.7))
        self.assertNoCosmic(wave, flux, (8344.1, 8344.8))
        
        #-- Test that absorption lines are still there
        self.assertAbsorptionLine(wave, flux, (5892.5, 5893.0), (5893.4, 5894.8), depth=0.6)
        self.assertAbsorptionLine(wave, flux, (7326.0, 7326.3), (7324.8, 7325.8), depth=0.5)
        
        #-- Test that the shift is correct
        f1 = flux[(wave>5883.71) & (wave<5883.94)]
        f2 = flux[(wave>5883.32) & (wave<5883.68)]
        self.assertTrue(np.min(f1) < np.min(f2), msg='Spectra NOT correctly shifted')
        f1 = flux[(wave>5914.00) & (wave<5914.30)]
        f2 = flux[(wave>5913.69) & (wave<5914.00)]
        self.assertTrue(np.min(f1) < np.min(f2), msg='Spectra NOT correctly shifted')
        
        
    #def testExt(self):
        #"""Merge original 51-orders hermes spectra"""
        #post_ext = '_HRF_OBJ_ext_CosmicsRemoved.fits'
        
class TestCase2SetHermesConfig(HermesTestCase):
    
    @classmethod
    def setUpClass(cls):
        #setup a directory for the test file if nessessary
        testdir = '/home/%s/temp/'%(getpass.getuser())
        if not os.path.exists(testdir):
            cls.dirExisted = False
            os.makedirs(testdir)
        else:
            cls.dirExisted = True
            
        cls.testdir = testdir
        cls.testfile = testdir+'hermesConfig.xml'
    
    def setUp(self):
        # create test hermesConfig.xml file
        config_file = "<hermes>\n"\
                      +"<Nights>/STER/mercator/hermes</Nights>\n"\
                      +"<CurrentNight>/STER/mercator/hermes/20110501</CurrentNight>\n"\
                      +"<Reduced>/STER/mercator/hermes/20110501/reduced</Reduced>\n"\
                      +"<AnalysesResults>/home/jorisv/hermesAnalyses</AnalysesResults>\n"\
                      +"<DebugPath>/home/jorisv/hermesDebug</DebugPath>\n"\
                      +"<LogFileName>/home/jorisv/hermesDebug/hermes.log</LogFileName>\n"\
                      +r"<LogFileFormater>%%(asctime)s  ::  %%(levelname)s  :: "\
                      +"%%(message)s</LogFileFormater>\n</hermes>"
        testfile = open(self.testfile, 'w')
        testfile.write(config_file)
        testfile.close()
        
        config = dict(Nights = '/STER/mercator/hermes',
                      CurrentNight = '/STER/mercator/hermes/20110501',
                      Reduced = '/STER/mercator/hermes/20110501/reduced',
                      AnalysesResults = '/home/jorisv/hermesAnalyses',
                      DebugPath = '/home/jorisv/hermesDebug',
                      LogFileName = '/home/jorisv/hermesDebug/hermes.log',
                      LogFileFormater = r"%%(asctime)s  ::  %%(levelname)s  :: %%(message)s")
        
        self.config = config
    
    def test1no_arguments(self):
        """catalogs.hermes write_hermesConfig no arguments"""
        
        old = hermes.write_hermesConfig(filename=self.testfile)
        
        msg = 'write_hermesConfig did not return the correct content of hermesConfig.xml'
        for key in old.keys():
            self.assertEqual(old[key], self.config[key], msg=msg)
        
        
    def test2correct_update(self):
        """catalogs.hermes write_hermesConfig update hermesConfig"""
        new_config = {'Nights' : '/STER/jorisv/hermesAnalysis',
                      'CurrentNight' : '/STER/jorisv/hermesAnalysis/BD+34.1543'}
        
        old = hermes.write_hermesConfig(filename=self.testfile, **new_config)
        new = hermes.write_hermesConfig(filename=self.testfile)
        
        msg = 'The new hermesConfig is not correctly written!'
        self.config.update(new_config)
        for key in new.keys():
            self.assertEqual(new[key], self.config[key], msg=msg) 
        
        
    def test3restore_old_config(self):
        """catalogs.hermes write_hermesConfig restore hermesConfig"""
        new_config = {'Nights' : '/STER/jorisv/hermesAnalysis',
                      'CurrentNight' : '/STER/jorisv/hermesAnalysis/BD+34.1543'}
        
        old = hermes.write_hermesConfig(filename=self.testfile, **new_config)
        new = hermes.write_hermesConfig(filename=self.testfile, **old)
        restored = hermes.write_hermesConfig(filename=self.testfile)
        
        msg = 'The new hermesConfig is not correctly restored!'
        for key in new.keys():
            self.assertEqual(restored[key], self.config[key], msg=msg) 
        
        
    @classmethod
    def tearDownClass(cls):
        #remove the test hermesConfig file
        os.remove(cls.testfile)
        
        #remove the test directory if nessessary
        if not cls.dirExisted:
            try:
                os.rmdir(cls.testdir)
            except Exception:
                print 'Error in removing directory'


class TestCase3ReadHermesVRVelocities(HermesTestCase):
    
    @classmethod
    def setUpClass(cls):
        #setup a directory for the test file if nessessary
        testdir = '/home/%s/temp/'%(getpass.getuser())
        if not os.path.exists(testdir):
            cls.dirExisted = False
            os.makedirs(testdir)
        else:
            cls.dirExisted = True
            
        cls.testdir = testdir
        cls.testfile = testdir + 'velocities.data'
        
        shutil.copy(cls.test_data_dir + 'catalogs_hermes_velocities.data',
                    cls.testfile)
    
    def setUp(self):
        # create test velocities.data file
        
        self.unseq = [2455351,2455425,2455620,2455639,2455649,2455665,2455937,2455351,
                      2455425,2455620,2455617,2455649,2455659]
        self.object = ['BD+29.3070', 'EC11031-1348', 'J1736+2806']
        
    def test1read_everything(self):
        """catalogs.hermes read_hermesVR_velocities basic read"""
        
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=False)
        
        msg = 'Not all data read from file, expected 13 lines, got {:}'.format(len(data))
        self.assertEqual(len(data), 13, msg=msg)

        msg = 'Not all expected unseq are present in data, expected:' \
              +'\n{:}\ngot:\n{:}'.format(set(self.unseq), set(data['unseq']))
        for sq in self.unseq: 
            self.assertTrue(sq in data['unseq'], msg=msg)
            
        msg = 'Not all expected objects are present, expected:' \
              +'\n{:}\ngot:\n{:}'.format(self.object, set(data['object']))
        for obj in self.object:
            self.assertTrue(obj in data['object'], msg=msg)
        
    def test2unique(self):
        """catalogs.hermes read_hermesVR_velocities unique output"""
        
        # return the earliest added lines
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=True,
                                               return_latest=False)
        
        msg = 'Amount of data not correct, expected 9 lines, got {:}'.format(len(data))
        self.assertEqual(len(data), 9, msg=msg)
        
        msg = 'All returned sequence numbers should be unique!, got: {:}'.format(data['unseq'])
        for sq in set(self.unseq):
            self.assertTrue(len(np.where(data['unseq'] == sq)) == 1, msg=msg)
        
        msg = 'When return_latest == False, Should return the first added lines, got the lastest!'
        self.assertEqual(data['vrad'][data['unseq'] == 2455351][0], -64.6618, msg=msg)
        
        # return the lastest added lines
        msg = 'When return_latest == True, Should return the latest added lines, got the first!'
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=True,
                                               return_latest=True)
        self.assertEqual(data['vrad'][data['unseq'] == 2455351][0], -66.5509, msg=msg)
        
        
    def test3list_unseq(self):
        """catalogs.hermes read_hermesVR_velocities list unseq"""
        
        unseq = [2455351, 2455425, 2455620]
        
        # return all
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=False,
                                               unseq=unseq)
        msg = 'Should return 2 lines per unseq provided\nGot:\t\t{:}\nExpected:\t{:}'.format( \
                data['unseq'], unseq+unseq)
        self.assertEqual(len(data['unseq']), 6, msg=msg)
        for sq in unseq:
            self.assertTrue(sq in data['unseq'], msg=msg)
            self.assertEqual(len(np.where(data['unseq']==sq)[0]), 2, msg=msg)
        
        # return unique
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=True,
                                               unseq=unseq)
        msg = 'Should return 1 lines per unseq provided\nGot:\t\t{:}\nExpected:\t{:}'.format( \
                data['unseq'], unseq)
        self.assertEqual(len(data['unseq']), 3, msg=msg)
        for sq in unseq:
            self.assertTrue(sq in data['unseq'], msg=msg)
            self.assertEqual(len(np.where(data['unseq']==sq)[0]), 1, msg=msg)
            
        # check returned order
        unseq = [2455937, 2455620, 2455425]
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=True,
                                               unseq=unseq)
        msg = 'The returned order should be the same as the provided list' \
              +'\nGot:\t\t{:}\nExpected:\t{:}'.format(data['unseq'], unseq)
        for i, sq in enumerate(unseq):
            self.assertEqual(data['unseq'][i], sq, msg=msg)
            
    def test4list_object(self):
        """catalogs.hermes read_hermesVR_velocities list object"""
        
        object = ['J1736+2806', 'EC11031-1348']
        unseq = [2455617, 2455649, 2455659, 2455665, 2455937]
        
        # return all objects
        data = hermes.read_hermesVR_velocities(filename=self.testfile, unique=False,
                                               object=object)
        msg = 'Should return 2 different object\nGot:\t\t{:}\nExpected:\t{:}'.format( \
                set(data['object']), object)
        self.assertEqual(len(set(data['object'])), 2, msg=msg)
        for ob in object:
            self.assertTrue(ob in data['object'], msg=msg)
        self.assertFalse('BD+29.3070' in data['object'], msg=msg)
        
        # return all lines belonging to those objects
        msg = 'Should return 5 lines in total\nGot:\t\t{:}\nExpected:\t{:}'.format( \
                set(data['unseq']), unseq)
        self.assertEqual(len(data['unseq']), 5, msg=msg)
        for sq in unseq:
            self.assertTrue(sq in data['unseq'], msg=msg)
        
        # return in the correct order
        object = ['J1736+2806','J1736+2806','J1736+2806','EC11031-1348','EC11031-1348']
        msg = 'Should return the objects in the correct order' \
              +'\nGot:\t\t{:}\nExpected:\t{:}'.format(data['object'], object)
        for i, ob in enumerate(object):
            self.assertEqual(data['object'][i], ob, msg=msg)
        
    @classmethod
    def tearDownClass(cls):
        #remove the test hermesConfig file
        os.remove(cls.testfile)
        
        #remove the test directory if nessessary
        if not cls.dirExisted:
            try:
                os.rmdir(cls.testdir)
            except Exception:
                print 'Error in removing directory'
        
        
class TestCase4ReadHermesVRAllCCF(HermesTestCase):
    """
    Tests regarding the function read_hermesVR_AllCCF of catalogs.hermes
    """
    
    @classmethod
    def setUpClass(cls):
        #setup a directory for the test file if nessessary
        testdir = '/home/%s/temp/'%(getpass.getuser())
        if not os.path.exists(testdir):
            cls.dirExisted = False
            os.makedirs(testdir)
        else:
            cls.dirExisted = True
        
        #create fake hermesAnalyses directory
        if not os.path.exists(testdir + 'hermesAnalyses/'):
            os.makedirs(testdir + 'hermesAnalyses/')
        
        #copy the allccf file
        shutil.copy(cls.test_data_dir + 'catalogs_hermes_00457640_AllCCF.fits',
                    testdir + 'hermesAnalyses/00457640_AllCCF.fits')
        
        cls.testdir = testdir
        cls.hermestestdir = testdir + 'hermesAnalyses/'
        cls.testfile = testdir + 'hermesAnalyses/00457640_AllCCF.fits'
        cls.unseq = 457640
    
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test1read_filename(self):
        """catalogs.hermes read_hermesVR_AllCCF read using filename"""
        
        hermes.HermesCCF.__init__ =mock.Mock(return_value=None) 
        
        ccf = hermes.read_hermesVR_AllCCF(unseq = self.testfile)
        ccf.__init__.assert_called_with(filename=self.testfile)
    
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test2read_unseq(self):
        """catalogs.hermes read_hermesVR_AllCCF read using unseq"""
        
        hermes.HermesCCF.__init__ =mock.Mock(return_value=None)
        hermes.hermesDir = self.testdir
        
        ccf = hermes.read_hermesVR_AllCCF(unseq = self.unseq)
        ccf.__init__.assert_called_with(filename=self.testfile)
            
    @classmethod
    def tearDownClass(cls):
        #remove the test hermesConfig file
        os.remove(cls.testfile)
        
        #remove the hermes testdir
        os.rmdir(cls.hermestestdir)
        
        #remove the test directory if nessessary
        if not cls.dirExisted:
            try:
                os.rmdir(cls.testdir)
            except Exception:
                print 'Error in removing directory'
                
        reload(hermes)
        
class TestCase5HermesCCF(HermesTestCase):
    """
    Tests regarding the class HermesCCF in catalogs.hermes
    """
    
    @classmethod
    def setUpClass(cls):
        #setup a directory for the test file if nessessary
        cls.testfile = cls.test_data_dir + 'catalogs_hermes_00457640_AllCCF.fits'
        cls.unseq = 457640
    
        cls.ccfobj = hermes.HermesCCF(filename=cls.testfile)
        
        print cls.ccfobj
        
    def test1init(self):
        """catalogs.hermes HermesCCF init"""
        
        ccf = self.ccfobj
        
        msg = 'Filename not correctly stored' \
              +"\nGot:\t\t{:}\nExpected:\t{:}".format( ccf.filename, self.testfile)
        self.assertEqual(ccf.filename, self.testfile)
        
        msg = "ccf object doesn't have all attributes"
        self.assertTrue(hasattr(ccf, 'vr'), msg=msg)
        self.assertTrue(hasattr(ccf, 'ccfs'), msg=msg)
        self.assertTrue(hasattr(ccf, 'groups'), msg=msg)
        
        msg = 'Header is not read correctly'
        self.assertEqual(ccf.header['object'], 'BD+42.3250', msg=msg)
        self.assertEqual(ccf.header['unseq'], 457640, msg=msg)
        
        msg = "VR list is not read correctly, got: len={:}, min={:}, max={:}".format(len(ccf.vr),
               ccf.vr[0], ccf.vr[-1])
        self.assertEqual(len(ccf.vr), 500, msg=msg)
        self.assertAlmostEqual(ccf.vr[0], -88.440071, places=3, msg=msg)
        self.assertAlmostEqual(ccf.vr[-1], 52.12719, places=3, msg=msg)
        
        msg = "Groups not correctly read, got:\n{:}".format(ccf.groups)
        names = ccf.groups['NAMES']
        self.assertEqual(names[0], '84-94 (viol)', msg=msg)
        self.assertEqual(names[3], '54-63       ', msg=msg)
        
        
    def test2ccf(self):
        """catalogs.hermes HermesCCF ccf()"""
        
        ccf = self.ccfobj
        
        # valid order
        msg = "ccf(80) should return valid ccf function"
        exp = [548.75676331, 567.31354311, 592.54494387, 589.58018404]
        got = ccf.ccf(80)[([10, 182, 362, 475],)]
        self.assertArrayAlmostEqual(exp, got, msg=msg, places=5)
        
        # invalid order
        msg = "invalid orders (<40 or >94) should raise IndexError"
        self.assertRaises(IndexError, ccf.ccf, 39)
        self.assertRaises(IndexError, ccf.ccf, 95)
        
    
    def test3combine_ccf(self):
        """catalogs.hermes HermesCCF combine_ccf()"""
        
        ccf = self.ccfobj
        
        # valid orders list
        msg = "combine_ccf([75,76,77,80,81,94]) should return valid ccf function"
        exp = [1811.14160156, 1871.66992188, 1902.21154785, 1985.3449707]
        got = ccf.combine_ccf([75,76,77,80,81,94])[([10, 182, 362, 475],)]
        self.assertArrayAlmostEqual(exp, got, msg=msg, places=5)
        
        # valid orders string
        msg = "combine_ccf('75:77,80,81,94') should return valid ccf function"
        exp = [1811.14160156, 1871.66992188, 1902.21154785, 1985.3449707]
        got = ccf.combine_ccf('75:77,80,81,94')[([10, 182, 362, 475],)]
        self.assertArrayAlmostEqual(exp, got, msg=msg, places=5)
        
        # invalid order list
        msg = "invalid orders (<40 or >94) should raise IndexError"
        self.assertRaises(IndexError, ccf.combine_ccf, [39,45,75,76])
        self.assertRaises(IndexError, ccf.combine_ccf, [45,75,76,95])
        
        # invalid order string
        msg = "invalid orders (<40 or >94) should raise IndexError"
        self.assertRaises(IndexError, ccf.combine_ccf, '39:45,75,76')
        self.assertRaises(IndexError, ccf.combine_ccf, '45:75,80:95')

class TestCase6RunHermesVR(HermesTestCase):
    """
    Integration Tests regarding run_hermesVR in catalogs.hermes
    """
    
    @classmethod
    def setUpClass(self):
        """ Mock all methods used by run_hermesVR """
        self.testfile = '/home/jorisv/IVS_python/ivs_test_data/catalogs_hermes_spectrum.fits'
        self.wvl_file = '/STER/mercator/hermes/20091130/reduced/00262831'
        
        hermes._subprocess_execute = mock.Mock(return_value=('', None, 0))
        hermes.write_hermesConfig = mock.Mock(return_value=dict(AnalysesResults='/home/'))
        os.remove = mock.Mock(return_value=None)
        os.mkdir = mock.Mock(return_value=None)
        os.rmdir = mock.Mock(return_value=None)
        shutil.copyfile = mock.Mock(return_value=None)
    
    def setUp(self):
        """ Create temporary directory """
        hermes.tempDir = self.make_testDir()
        
    def tearDown(self):
        """ Delete temporary directory """
        self.delete_testDir()
    
    @classmethod
    def tearDownClass(self):
        """ restore all mocked methods """
        reload(os)
        reload(shutil)
        reload(hermes)
        
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test1create_delete_dir(self):
        """catalogs.hermes run_hermesVR creating/deleting temp directory"""
        
        
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='release',
                                                     unseq=500)
        
        #-- create temp directories
        os.mkdir.assert_any_call(hermes.tempDir+'hermesvr/')
        os.mkdir.assert_any_call(hermes.tempDir+'hermesvr/reduced/')
        
        
        
        #-- delete temp directories
        os.rmdir.assert_any_call(hermes.tempDir+'hermesvr/')
        os.rmdir.assert_any_call(hermes.tempDir+'hermesvr/reduced/')        
        
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test2copy_input(self):
        """catalogs.hermes run_hermesVR copy input files"""
        
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='release',
                                                     unseq=500)
        
        #-- copy the input file
        shutil.copyfile.assert_called_with(self.testfile, 
                                           hermes.tempDir + 'hermesvr/reduced/500_HRF_OBJ_ext.fits')
        
        #-- delete the file again
        os.remove.assert_called_with(hermes.tempDir+'hermesvr/reduced/500_HRF_OBJ_ext.fits')
    
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test3hermes_config(self):
        """catalogs.hermes run_hermesVR write correct hermesConfig"""
        
        #-- release version
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='release',
                                                     unseq=500)
        
        hermes.write_hermesConfig.assert_any_call(Nights = hermes.tempDir, 
                                            CurrentNight = hermes.tempDir + 'hermesvr/',
                                            Reduced = hermes.tempDir)
        
        #-- trunk version
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='trunk',
                                                     unseq=500)
        
        hermes.write_hermesConfig.assert_any_call(Nights = hermes.tempDir + 'hermesvr/', 
                                            CurrentNight = hermes.tempDir + 'hermesvr/',
                                            Reduced = hermes.tempDir)
    
    @unittest.skipIf(noMock, "Mock not installed")
    def test2defaults(self):
        """catalogs.hermes run_hermesVR correct hermesVR command"""
        
        #-- release version
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='release')
        
        cmd = ['python', '/STER/mercator/mercator/Hermes/releases/hermes5/pipeline/run/hermesVR.py',
               '-i', '262834', '-w', '/STER/mercator/hermes/20091130/reduced/00262831']
        hermes._subprocess_execute.assert_called_with(cmd, 500)
        
        #-- trunk version
        unseq, out, returncode = hermes.run_hermesVR(self.testfile, wvl_file=self.wvl_file, version='trunk',
                                                     timeout=150)
        
        cmd = ['python', '/STER/mercator/mercator/Hermes/trunk/hermes/pipeline/run/hermesVR.py',
               '-i', '262834', '-w', '/STER/mercator/hermes/20091130/reduced/00262831', '-f']
        hermes._subprocess_execute.assert_called_with(cmd, 150)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
        
