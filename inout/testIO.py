import os
import h5py
import numpy as np
from ivs.inout import hdf5

import unittest

class HDF5TestCase(unittest.TestCase):

    def testWriteDict(self):
        """ inout.hdf5.write_dict() """
        data = {}
        data['number'] = 12.30
        data['list'] = [45, 78.36]
        data['string'] = 'This is a piece of text'

        set1 = {}
        set1['grid'] = np.random.uniform(size=(10,10))
        set1['x'] =  np.linspace(0,10,num=10)
        set1['y'] = np.linspace(-5,5,num=10)
        set1['subset'] = {'w1':10.0, 'w2':12.0, 'w3':13.0}
        data['set1'] = set1

        set2 = {}
        set2['grid'] = np.random.uniform(size=(10,10))
        set2['x'] =  np.linspace(-10,0,num=10)
        set2['y'] = np.linspace(0,5,num=10)
        set2['subset'] = {'w1':1, 'w2':2, 'w3':3, 'label':'again some text'}
        data['set2'] = set2

        hdf5.write_dict(data, 'test.hdf5', update=False, attr_types=[float, int, str, list])

        hdf = h5py.File('test.hdf5', 'r')

        self.assertTrue('set1' in hdf)
        self.assertTrue('set2' in hdf)
        self.assertTrue('number' in hdf.attrs)
        self.assertTrue('string' in hdf.attrs)
        self.assertTrue('list' in hdf.attrs)

        self.assertEqual(hdf.attrs['number'], data['number'])
        self.assertEqual(hdf.attrs['string'], data['string'])
        self.assertListEqual(hdf.attrs['list'].tolist(), data['list'])

        self.assertTrue(np.all(hdf['set1']['grid'].value == set1['grid']))
        self.assertTrue(np.all(hdf['set1']['x'].value == set1['x']))
        self.assertTrue(np.all(hdf['set1']['y'].value == set1['y']))

        self.assertEqual(hdf['set1']['subset'].attrs['w1'], 10.0)
        self.assertEqual(hdf['set1']['subset'].attrs['w2'], 12.0)
        self.assertEqual(hdf['set1']['subset'].attrs['w3'], 13.0)

        self.assertTrue(np.all(hdf['set2']['grid'].value == set2['grid']))
        self.assertTrue(np.all(hdf['set2']['x'].value == set2['x']))
        self.assertTrue(np.all(hdf['set2']['y'].value == set2['y']))

        self.assertEqual(hdf['set2']['subset'].attrs['w1'], 1)
        self.assertEqual(hdf['set2']['subset'].attrs['w2'], 2)
        self.assertEqual(hdf['set2']['subset'].attrs['w3'], 3)

        hdf.close()

        if os.path.isfile('test.hdf5'):
            os.remove('test.hdf5')

    def testReadDict(self):
        """ inout.hdf5.read2dict() """

        if os.path.isfile('test.hdf5'):
            os.remove('test.hdf5')
        hdf = h5py.File('test.hdf5', 'w')

        hdf.create_group('set1')
        hdf.create_group('set2')
        hdf.create_dataset('list', data=[12,124.54])

        grid1 = np.random.uniform(size=(10,10))
        hdf['set1'].create_dataset('grid', data = grid1)
        hdf['set1'].attrs['w1'] = 'this is a label'
        hdf['set1'].create_group('subset')
        hdf['set1']['subset'].attrs['w1'] = 15.2
        hdf['set1']['subset'].attrs['w2'] = -1000000000000
        hdf['set1']['subset'].attrs['w3'] = 'yet another text string'

        grid2 = np.random.uniform(size=(10,10))
        hdf['set2'].create_dataset('grid', data = grid2)

        hdf.close()

        data = hdf5.read2dict('test.hdf5')

        self.assertTrue('set1' in data)
        self.assertTrue('set2' in data)
        self.assertTrue('list' in data)
        self.assertListEqual(data['list'].tolist(), [12,124.54])

        self.assertTrue(np.all(data['set1']['grid'] == grid1))
        self.assertEqual(data['set1']['w1'], 'this is a label')

        self.assertEqual(data['set1']['subset']['w1'], 15.2)
        self.assertEqual(data['set1']['subset']['w2'], -1000000000000)
        self.assertEqual(data['set1']['subset']['w3'], 'yet another text string')
        self.assertTrue(np.all(data['set2']['grid'] == grid2))

        if os.path.isfile('test.hdf5'):
            os.remove('test.hdf5')



