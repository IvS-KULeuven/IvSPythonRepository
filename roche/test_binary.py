"""Module test_binary - """

import binary
import logging

logging.basicConfig(level=logging.DEBUG)
no_assert = 'There is currently no Assert check implemented'

# global test values
x = 1
y = 1
z = 1
r = 1
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

            log.debug('For case = ' + str(norm))
            log.debug(no_assert)
            # log.debug()

    def test_binary_roche_surface_gravity(self):
        log = logging.getLogger('roche_surface_gravity')

        for norm in booleans:
            binary.binary_roche_surface_gravity(x, y, z, d, omega, M1, M2,
                                                norm)
            log.debug('For case = ' + str(norm))
            log.debug(no_assert)

    def test_get_binary_roche_radius(self):
        log = logging.getLogger('get_binary_roche_radius')

        binary.get_binary_roche_radius(theta, phi, Phi, q, d, F, r_pole=None)
        log.debug(no_assert)
