# -*- coding: utf-8 -*-
"""
Classe SSE

@author: mathe
"""
import sys
import numpy as np
from result import Result 

class SSE:
    """ 
    Classe para montar problema de calculo do SSE
    """
    
    def __init__(self,material,ensaios,experimentos):
        """
        Construtor classe SSE
        """
        n_ensaios=len(ensaios)
        n_experimentos=len(experimentos)
        if n_ensaios != n_experimentos:
            sys.exit('Fatal Error: n. ensaios e experimentos não conferem')
        
        self.material=material
        self.ensaios=ensaios
        self.experimentos=experimentos
        self.sse0=1
        self.resultados=[None]*n_ensaios
        
        for i in range(n_ensaios):
            self.resultados[i]=Result()
            self.resultados[i].set_Result(self.material,self.ensaios[i])
            self.ensaios[i].criar_mapeamento_experimentos(self.experimentos[i].t)
        
#-------------------------------------------------------------#
    
    def set_normal(self,sse0):
        self.sse0=sse0

#-------------------------------------------------------------#
    
    def avaliar(self,x):
        """
        Metodo para avaliação do SSE 
        """
        self.material.set_prop(x)
        n_ensaios=len(self.ensaios)
        self.sse_ensaio=np.zeros(n_ensaios)
        
        for i in range(n_ensaios):
            self.resultados[i],fail=self.ensaios[i].run(self.material,self.resultados[i])
        
            if fail:
                sys.exit('Fatal error: Erro do Ensaio')
        
            self.sse_ensaio[i]=self.ensaios[i].calc_sse(self.experimentos[i],self.resultados[i])
        
        sse=np.sum(self.sse_ensaio)/self.sse0
        return sse
        
#-------------------------------------------------------------#    
            