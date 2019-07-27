#! /usr/bin/env python2.6
# -*- coding: utf-8 -*-

import integration
import math

integrator = integration.Integration(0,math.pi,100)

print 'trapezoidal', integrator.trapezoidal_rule(integrator.test_functional,1,0,0)
print 'rectangle', integrator.rectangle_rule(integrator.test_functional,1,0,0)
print 'basic', integrator.riemman_rule(integrator.test_functional,1,0,0)

print integrator.trapezoidal_rule(integrator.test_functional,1,0,0) - integrator.rectangle_rule(integrator.test_functional,1,0,0)
