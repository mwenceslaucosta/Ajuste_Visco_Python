"""
Teste modelo em tracao
Comparado com resultados obtidos pelo Matlab
Utilizar rotina do matlab: Teste_Zener3D_classico_com_implementaca_python
"""
import math 
import numpy as np
from matplotlib import pyplot as plt
from modulo_ensaios import Ensaio_tracao
from modulo_modelosZener3D_infinitesimal import ModeloZener3D_classico_infini
from result import Result


# Instanceando classe ensaio de tracao
ensaio=Ensaio_tracao()

# Definindo valores dos atributos 
tensile={}
tensile['def_rate']=1.0E-2
tensile['def_max']=0.1
tensile['step_tempo']=300
ensaio.set_ensaio(tensile)

# obtendo tempo maximo para ensaio numerico 
max_n_time=ensaio.get_n_ensaio()
config={}
config['n_time']=max_n_time

# Instanceando classe modelo material 
config['k_bracos']=2
config['poisson']=0.46
material=ModeloZener3D_classico_infini(config)

# Parametos materiais de teste 
x=np.array([3.84, 38.46, math.log10(1), 38.46/2, math.log10(0.5)])

resultados=Result()
resultados.set_Result(material,ensaio)
material.set_prop(x)
resultados,fail=ensaio.run(material,resultados)

# Plotando resultados 
fig1=plt.figure()
plt.plot(resultados.e[0,:],resultados.stress[0,:])