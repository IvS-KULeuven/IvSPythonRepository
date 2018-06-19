"""Module test_binary - """

import binary
import logging
import numpy as np

logging.basicConfig(level=logging.DEBUG)
no_assert = 'There is currently no Assert check implemented'

# global test values
x = 1
y = 1
z = 1
r = 1
r_pole = 1
theta = 0
phi = 0
Phi = 0
d = 1
F = 1
omega = 1
M1 = 1
M2 = 1
q = 1
booleans = [True, False]
star = np.rec.fromarrays([[1,1,1,1]]*16, names=['theta', 'phi', 'r', 'x',
                                                  'y', 'z', 'vx', 'vy', 'vz',
                                                  'gravx', 'gravy', 'gravz',
                                                  'grav', 'areas', 'teff',
                                                  'flux']
                         )


class TestRochePotential(object):

    def test_print_tester(self):
        name = 'abc'
        assert binary.print_tester(name) == 'ABC'

    def test_binary_roche_potential(self):
        log = logging.getLogger('roche_potential')

        binary.binary_roche_potential(r, theta, phi, Phi, q, d, F)
        log.debug(no_assert)
        # print(no_assert)

    def test_binary_roche_potential_gradient(self):
        log = logging.getLogger('roche_potential_gradient')

        for norm in booleans:
            binary.binary_roche_potential_gradient(x, y, z, q, d, F, norm)

            log.debug('For case norm = {}'.format(norm))
            log.debug(no_assert)
            # log.debug()

    def test_binary_roche_surface_gravity(self):
        log = logging.getLogger('roche_surface_gravity')

        for norm in booleans:
            binary.binary_roche_surface_gravity(x, y, z, d, omega, M1, M2,
                                                norm)
            log.debug('For case norm = {}'.format(norm))
            log.debug(no_assert)

    def test_get_binary_roche_radius(self):
        log = logging.getLogger('get_binary_roche_radius')

        binary.get_binary_roche_radius(theta, phi, Phi, q, d, F, r_pole)
        log.debug(no_assert)

    def test_reflection_effect(self):
        log = logging.getLogger('reflection_effect')

        binary.reflection_effect(star, star, theta, phi)
        lof.debug(no_assert)

    def test_spectral_synthesis(test):
        log = logging.getLogger('spectral_synthesis')

        binary.spectral_synthesis(star)
        lof.debug(no_assert)

    def test_binary_light_curve_synthesis(test):
        log = logging.getLogger('binary_light_curve_synthesis')
        parameters = {}
        parameters['Tpole1'] = 1     # Primary Polar temperature   [K]
        parameters['Tpole2'] = 1    # Secondary Polar temperature [K]
        parameters['P'] = 1             # Period [days]
        parameters['asini'] = 1        # total semi-major axis*sini  [AU]
        parameters['Phi1'] = 0          # Gravitational potential of primary [-]
        parameters['Phi2'] = 0

        binary.binary_light_curve_synthesis(parameters=parameters)
        lof.debug(no_assert)
