# -*- coding: utf-8 -*-
"""
Main ajuste de parâmetros dados em fluência-recuperação Brant 
Utiliza ModeloZener3D_classico_infini

@author: mathe
"""
import numpy as np 
import matplotlib.pyplot as plt
import time
from modulo_modelosZener3D_infinitesimal import ModeloZener3D_classico_infini
from modulo_dados_experimentais import Brant_creep_reco_UHMWPE
from protocolo_otimizacao import Protocolo_De_Otimizacao
from modulo_sse import SSE 

config={}

#Modelo Material
config['k_bracos']=1
config['poisson']=0.46
material=ModeloZener3D_classico_infini(config)

delta_t=100

#Instanceando_classe_dados_experimentais utilizada no problema de ajuste
vetor_ensaios=np.array([1,2,3])
creep_reco=Brant_creep_reco_UHMWPE()
ensaios,experimentos=creep_reco.set_ensaios_experimentos(vetor_ensaios,delta_t)

#Definindo tempo maximo numerico
max_n_time=0
for i in range(len(ensaios)):
    max_n_time=np.max([max_n_time,ensaios[i].get_n_ensaio()])

n_time=max_n_time
material.value_n_time(n_time)

#Configuração problema de otimização
#Restrições 
config_PSO={}
config_PSO['LB']=np.array([0,1,-4,1,-4])
config_PSO['UB']=np.array([4.3, 4.3, 6, 4.3, 6])
x0=np.array([2.0130,2.1575,3.6045])  
#Chamada Protocolo de otimização
sse=SSE(material,ensaios,experimentos) 
t = time.time()
for i in range(5):
    print(i)
    a=sse.avaliar(x0)
elapsed = time.time() - t
#protocolo=Protocolo_De_Otimizacao(config_PSO)
#result=protocolo.run_protocolo(material,ensaios,experimentos)

