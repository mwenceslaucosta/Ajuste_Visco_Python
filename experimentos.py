# -*- coding: utf-8 -*-
"""
Modulo para importar dados experimentais

@author: mathe
"""
import numpy as np 
import sys
class Experimentos: 
    """ 
    classe Experimentos 
    """
    def set_experimento(self,file_name):
        tabela=np.loadtxt(file_name,delimiter=",")
        number_cols=tabela.shape[1]
        if (number_cols < 2):
            sys.exit('Falal error: Número de colunas não conferem')  
        
        self.t=tabela[:,0]
        self.e=tabela[:,1]
        
        if (number_cols>2):
            self.stress=tabela[:,2]
        
        
