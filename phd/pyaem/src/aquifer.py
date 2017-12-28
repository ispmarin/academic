'''
Created on 07/03/2011

@author: ispmarin
'''

class AquiferData():
    '''
    classdocs
    '''


    def __init__(self,sbase,sheight,sk_ext,sreference_head,sreference_point):
        '''
        Constructor
        '''
        self.base = sbase
        self.height = sheight
        self.k_ext = sk_ext
        self.reference_head = sreference_head
        self.reference_point = sreference_point
        
        #print( 'Aquifer ', 'base', self.base,'height', self.height,'k_ext', self.k_ext,'reference head', self.reference_head,'reference_point', self.reference_point)
        

        
    