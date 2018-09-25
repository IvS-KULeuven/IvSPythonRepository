"""
Module to test builder.py

These units test are writting to be used by pytest.
"""

import builder


class TestSEDobject(object):

    def setup_method(self):
        ID = 'HD_35155'
        print('Setting up SED object with ID = {}'.format(ID))
        self.sed = builder.SED(ID)

    def test_init(self):
        assert isinstance(self.sed, builder.SED)
