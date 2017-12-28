'''
Created on 29/03/2011

@author: ispmarin
'''

import math
#import numpy as np
    

class CreateEllipse(object):
    '''
    classdocs
    '''

    def __init__(self, sa, sb, sangle, scenter, ssubdivisions):
        '''
        Constructor
        '''
        self.a = sa
        self.b = sb
        self.angle = sangle #* math.pi / 180.0
        self.center = scenter
        self.subdivisions = ssubdivisions
        self.vertices = []
        
        
    def __repr__(self):
        return "Ellipse(a=%s, b=%s, angle=%s, center=%s, subdivisions=%s)" % \
                                                  (self.a, 
                                                   self.b, 
                                                   self.angle,
                                                   self.center,
                                                   self.subdivisions)
        
    @classmethod
    def by_semi_axis(cls, sa, sb, sangle, scenter, ssubdivisions):
        return cls(sa, sb, sangle, scenter, ssubdivisions)
    
    @classmethod
    def by_box(cls, x_inf, x_sup, sb, ssubdivisions):
        
        a = 0.5 * math.sqrt(pow((x_sup.real) - (x_inf.real), 2) + \
        pow((x_sup.imag) - (x_inf.imag), 2))
        
        b = sb * 0.5 #NOTE: this may be wrong
        
        angle = math.atan2(-((x_sup.imag) - (x_inf.imag)), \
        -((x_sup.real) - (x_inf.real))) + math.pi
        center = complex((x_inf.real) + ((x_sup.real) - \
        (x_inf.real)) * 0.5, (x_inf.imag) + ((x_sup.imag) - (x_inf.imag)) * 0.5)
        
        return cls(a,b,angle,center,ssubdivisions)
    
    def print_ellipse(self):
        for vertice in self.vertices:
            print vertice.real, vertice.imag
        
    def create_ellipse(self):
        
        step = 2.0 * math.pi /self.subdivisions

        for i in xrange(0,self.subdivisions): 
            self.vertices.append( complex ((self.center.real) + self.a * math.cos(i*step) * math.cos(self.angle) - self.b
                    * math.sin(i*step) * math.sin(self.angle), (self.center.imag) + self.a * math.cos(i*step) * math.sin(
                    self.angle) + self.b * math.sin(i*step) * math.cos(self.angle)))
        
        #self.print_ellipse()

               
