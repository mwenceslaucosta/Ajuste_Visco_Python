"""
Teste modelo em fluencia-recuperacao
Comparado com resultados obtidos pelo Matlab
Utilizar rotina do matlab: Teste_Zener3D_classico_com_implementaca_python
"""
import math 
import numpy as np
from matplotlib import pyplot as plt
from modulo_ensaios import Ensaio_fluencia_recuperacao
from modulo_modelosZener3D_infinitesimal import ModeloZener3D_classico_infini
from result import Result
import time

# Instanceando classe ensaio de tracao
ensaio=Ensaio_fluencia_recuperacao()

# Definindo valores dos atributos 
fluencia={}
fluencia['stress']=4
fluencia['sigma_recovery']=0
fluencia['tf']=165E3
fluencia['delta_t']=100
fluencia['t_carregamento']=[3, 2*60, 24*60*60]
ensaio.set_ensaio(fluencia)

# obtendo tempo maximo para ensaio numerico 
max_n_time=ensaio.get_n_ensaio()
config={}
n_time=max_n_time

# Instanceando classe modelo material 
config['k_bracos']=2
config['poisson']=0.46
material=ModeloZener3D_classico_infini(config)

# Parametos materiais de teste 
x=np.array([0.959474268362036,1.51569752353631,5.18653621382354,1.97476513530482,3.05140823297537])

resultados=Result()
resultados.set_Result(material,ensaio)
material.set_prop(x)
t = time.time()
for i in range(30):
    print(i)
    resultados,fail=ensaio.run(material,resultados)
elapsed = time.time() - t

# Plotando resultados 
fig1=plt.figure()
plt.plot(resultados.t,resultados.e[0,:])