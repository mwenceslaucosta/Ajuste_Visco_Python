# -*- coding: utf-8 -*-
"""
Modulo com Resultados experimentais utilizados nos ajustes 

@author: mathe
"""
import numpy as np 
from modulo_ensaios import Ensaio_fluencia_recuperacao
from experimentos import Experimentos 

class Brant_creep_reco_UHMWPE:
    """
    Classe para importar dados experimentais em fluencia-recuperação Brant 
    e configurar ensaios numericos para estes dados.
    """
    def set_ensaios_experimentos(self,vetor_ensaios,delta_t):
        """
        Metodo para setar experimentos e ensaios numericos em fluencia-recuperacao
        Brant
        """
        n_ensaios=len(vetor_ensaios)
        ensaios=[None]*n_ensaios
        experimentos=[None]*n_ensaios
        fluencia={}
        for i in range(n_ensaios):
            if (vetor_ensaios[i]==1): 
                ensaios[i]=Ensaio_fluencia_recuperacao()
                experimentos[i]=Experimentos()
                fluencia['stress']=4
                fluencia['sigma_recovery']=0.07
                fluencia['tf']=165E3
                fluencia['delta_t']=delta_t
                fluencia['t_carregamento']=np.array([3, 2*60,24*60*60])
                ensaios[i].set_ensaio(fluencia)
                experimentos[i].set_experimento('Brant_creep_reco_4MPa.csv')

            if (vetor_ensaios[i]==2):
                ensaios[i]=Ensaio_fluencia_recuperacao()
                experimentos[i]=Experimentos()               
                fluencia['stress']=8
                fluencia['sigma_recovery']=0.07
                fluencia['tf']=165E3
                fluencia['delta_t']=delta_t
                fluencia['t_carregamento']=np.array([3, 4*60,24*60*60])
                ensaios[i].set_ensaio(fluencia)
                experimentos[i].set_experimento('Brant_creep_reco_8MPa.csv') 

            if (vetor_ensaios[i]==3):
                ensaios[i]=Ensaio_fluencia_recuperacao()
                experimentos[i]=Experimentos()               
                fluencia['stress']=16
                fluencia['sigma_recovery']=0.07
                fluencia['tf']=165E3
                fluencia['delta_t']=delta_t
                fluencia['t_carregamento']=np.array([3, 8*60,24*60*60])
                ensaios[i].set_ensaio(fluencia)
                experimentos[i].set_experimento('Brant_creep_reco_16MPa.csv')
        
        return ensaios,experimentos
                
                
        
        
        