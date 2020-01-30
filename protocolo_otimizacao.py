# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 19:54:14 2020

@author: mathe
"""
from modulo_sse import SSE
#import pyswarms as PSO
import numpy as np 
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import BFGS
from pyswarm import pso as PSO

class Protocolo_De_Otimizacao: 
    """ 
    Classe protocolo de otimizacao. Chama problema de otimizacao de parametros
    """
#-----------------------------------------------------------------------------#
    def __init__(self,config_PSO):
        """
        Construtor da classe Protocolo_De_Otimizacao
        """
        if config_PSO.get('SwarmSize') is None: 
            self.SwarmSize=150
        else: 
            self.SwarmSize=config_PSO['SwarmSize']
        if config_PSO.get('Inertia') is None:
            self.w=0.5
        else: 
            self.w=config_PSO.get['Inertia']
        if config_PSO.get('SelfAdjustment') is None: 
            self.c1=1.3
        else: 
            self.c1=config_PSO['SelfAdjustment']
        if config_PSO.get('SocialAdjustment') is None: 
            self.c2=1.3
        else: 
            self.c2=config_PSO['SocialAdjustment']
        if config_PSO.get('MaxIterations') is None: 
            self.MaxIter=5
        else: 
            self.MaxIter=config_PSO['MaxIterations']
        if config_PSO.get('n_PSO') is None: 
            self.n_PSO=1
        else: 
            self.n_PSO=config_PSO['n_PSO']
        
        self.max_bound=config_PSO['UB']
        self.min_bound=config_PSO['LB']
        self.gradient_tol=1E-10

#-----------------------------------------------------------------------------#
    def run_protocolo(self,material,ensaios,experimentos):
        """
        Método RUN para protocolo de otimização
        """
        #flag_particle - Utilizado para entrar com cada particula em separado 
        # na funcao objetivo, caso contrario a biblioteca PSO entraria na funçao
        # objetivo (avaliar.sse) com uma matriz contendo todos os vetores de projeto.
        #flag_particle=True resolve para cada particula (individuo) em separado
        
        options_PSO={'c1':self.c1,'c2':self.c2,'w':self.w}
        dim=len(self.max_bound)        
        bounds_PSO = (self.min_bound, self.max_bound)
        
        result={}
        result['SSE_PSO']=[None]*self.n_PSO
        cost_PSO=np.zeros(self.n_PSO)
        x_PSO=np.zeros((self.n_PSO,dim))
        
        for i in range(self.n_PSO):
            print('Iniciada ', i+1,'° iteração PSO')
            result['SSE_PSO'][i]=SSE(material,ensaios,experimentos)      
            #Instanciando PSO da biblioteca pyswarms
            #optimizer = PSO.single.GlobalBestPSO(n_particles=self.SwarmSize,\
            #            dimensions=dim, options=options_PSO, bounds=bounds_PSO)
            
            #Chamada efetiva otimizador PSO
                
            #flag_particle=True                                            
            #cost_PSO[i], x_PSO[i,:] = optimizer.optimize\
            #    (result['SSE_PSO'][i].avaliar, iters=self.MaxIter,flag_particle=flag_particle)
            
            x_PSO[i,:],cost_PSO[i]=PSO(result['SSE_PSO'][i].avaliar,lb=self.min_bound,\
                        ub=self.max_bound,swarmsize=self.SwarmSize, omega=self.w,\
                  phip=self.c1, phig=self.c2, maxiter=self.MaxIter,debug=True)
                
            #Salvando dados de tensao e deforação dos ensaios para cada PSO    
            tabela_strain=np.zeros((ensaios[1].n,len(ensaios)))
            tabela_stress=np.zeros((ensaios[1].n,len(ensaios)))
                        
            for j in range(len(ensaios)):
                tabela_strain[:,i] = result['SSE_PSO'][i].resultados[j].e[1,:]
                tabela_stress[:,i] = result['SSE_PSO'][i].resultados[j].stress[1,:]
            np.savetxt('stess_PSO' + str(i) + '.csv',tabela_stress,delimiter=",")
            np.savetxt('strain_PSO' + str(i) +'.csv',tabela_strain,delimiter=",")
        
        #Salvando vetor tempo 
        np.savetxt('time.csv',result['SSE_PSO'][1].resultados[1].t)
        result['cost_PSO'][:]=cost_PSO
        result['x_PSO'][:,:]=x_PSO
        #Definindo melhor iteração PSO
        result['index_min_PSO']=np.argmin(cost_PSO)
        result['min_SSE_PSO']=cost_PSO[index_min_PSO]
        
        
        x0=x_PSO[result['index_min_PSO'],:]
        """
        x0=np.array([0.969779497099175,1.51324756484957,5.18154144137473,1.97509730202257,3.04983344905473])
        LB=x0-1.5
        UB=x0+1.5
        
        if (material.tipo_modelo=='fracionario'):
            cont=0
            for i in range(material.k_bracos):
                LB[cont+3]=0.01
                LB[cont+3]=1
                cont+=3
        
        bounds = Bounds(LB,UB)  
        result['Hibrido']=SSE(material,ensaios,experimentos)  
        print('Hibrido iniciado')
        options_trust_reg={'disp': True,'verbose': 2}
        #options_trust_reg={'gtol':self.gradient_tol,'xtol':self.gradient_tol,'disp': True,'verbose': 2}
        res_hib = minimize(result['Hibrido'].avaliar, x0, method='trust-constr',\
                   hess=BFGS(),options=options_trust_reg, bounds=bounds)            
        #options={'ftol': 1E-10,'disp': True}
        #res = minimize(result['Hibrido'].avaliar, x0, method='SLSQP',
                  #options=options,bounds=bounds)
        result['x_otimo']=res_hib.x
        result['f_otimo']=res_hib.fun
        #Salvando dados de tensao e deforação dos ensaios para Hibrido
        
        tabela_strain_hibr=np.zeros((ensaios[1].n,len(ensaios)))
        tabela_stres_hibr=np.zeros((ensaios[1].n,len(ensaios)))
        for j in range(len(ensaios)):
            tabela_strain_hibr[:,i] = result['Hibrido'].resultados[j].e[1,:]
            tabela_stres_hibr[:,i] = result['Hibrido'].resultados[j].stress[1,:]
        np.savetxt('stess_Hibrido.csv',tabela_stress,delimiter=",")
        np.savetxt('strain_Hibrido.csv',tabela_strain,delimiter=",")

        #Salvando vetor de propriedades
        x_otimos=np.array([[x_PSO],[res_hib.x]])
        np.savetxt('X_otimos.csv',x_otimos,delimiter=",")
        """
        return result
    
    
    
    
    

#-----------------------------------------------------------------------------#    