# -*- coding: utf-8 -*-
import numpy as np

def Ensaio_ficticio_fluencia(config,delta_t):
    """
    Funcao para gerar informações necessarias para criar ensaio ficticio de fluencia
    """
    vetor_ensaios=config['vetor_ensaios']
    n_ensaios=len(vetor_ensaios)
    ensaios=[None]*n_ensaios
        
    for i in range(n_ensaios):
        # Caso 1 - Ensaio 1
        if vetor_ensaios[i]==1:
            fluencia['stress']=4
            fluencia['tf']=86E3
            fluencia['delta_t']=delta_t
            fluencia['t_carregamento']=np.([1, 2*60, 24*60*60])
            ensaios[i]=config['fun_ensaio']
            ensaios[i].set_ensaio(fluencia)
       # Caso 2 - Ensaio 2
        if vetor_ensaios[i]==2:
            fluencia['stress']=8
            fluencia['tf']=86E3
            fluencia['delta_t']=delta_t
            fluencia['t_carregamento']=np.([1, 4*60, 24*60*60])
            ensaios[i]=config['fun_ensaio']
            ensaios[i].set_ensaio(fluencia)
        # Caso 3 - Ensaio 3
        if vetor_ensaios[i]==3:
            fluencia['stress']=16
            fluencia['tf']=86E3
            fluencia['delta_t']=delta_t
            fluencia['t_carregamento']=np.([1, 8*60, 24*60*60])
            ensaios[i]=config['fun_ensaio']
            ensaios[i].set_ensaio(fluencia)
        
    return ensaios
