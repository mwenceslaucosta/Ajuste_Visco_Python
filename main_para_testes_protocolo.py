# -*- coding: utf-8 -*-
"""
Main para Testes do protocolo 
"""
from modulo_modelosZener3D_infinitesimal import ModeloZener3D_classico_infini
from modulo_ensaios import Ensaio_fluencia
config=[]
# Adicionar set_path

# Configuracao do problema para analise do intervalo de tempo 
config['fun_ensaio']=Ensaio_fluencia
config['k_bracos']=2
config['poisson']=0.46
config['fun_material']=ModeloZener3D_classico_infini(config)

# Restricoes do problema de otimizacao 
config['LB']=np.array([1, 1, -4, 1, -4])
config['UB']=np.array(4.3, 4.3, 6, 4.3, 6])

# Parametros de referencia para gerar dados "experimentais numericos"

# Fazer problema de verificacao delta_t
delta_t=100

# Configuracao problema de otimizacao 
config['vetor_ensaios']=np.
