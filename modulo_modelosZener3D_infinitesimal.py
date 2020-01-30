# -*- coding: utf-8 -*-
"""
Modulo com modelos Zener 3D clássicos

@author: mathe
"""

import numpy as np 
import math
from modelos import Modelos 

#-------------------------------------------------------------#     
#--------------------Zener 3D Classico------------------------#
#-------------------------------------------------------------#
class ModeloZener3D_classico_infini(Modelos):
    """
    Modelo Zener 3D classico em cinematica infinitesimal.
    """
    def __init__(self,config):
        """
        Construtor 
        """
        self.k_bracos=config['k_bracos']
        self.dimension=6
        self.poisson=config['poisson']
        self.tipo_modelo='classico'
#-------------------------------------------------------------     
    def value_n_time(self,n_time):
        """
        armazena vetor tempo 
        """
        self.n_time=n_time
#-------------------------------------------------------------         
    def get_size_var_inter(self):
        """
        Fornece numero de variaveis internas e dimensao
        """
        size=np.zeros(2,dtype=int)
        size[0]=self.dimension*self.k_bracos
        size[1]=self.dimension 
        return size 
#-------------------------------------------------------------     
    def set_prop(self,x):
        """
        Atualiza vetor de atributos do modelo 
        """
        self.G_inf=10**x[0]
        self.G=np.zeros(self.k_bracos)
        self.tau=np.zeros(self.k_bracos)
        cont_1=1
        
        self.G_0=self.G_inf
#-------------------------------------------------------------      
        for i in range(self.k_bracos):

            self.G[i]=10**x[cont_1+i]
            self.tau[i]=10**x[cont_1+1+i]
            self.G_0+=self.G[i]
            cont_1+=1
        
        self.K_inf=2*self.G_0*(1+self.poisson)/(3*(1-2*self.poisson))


#-------------------------------------------------------------         
    def mat_solve(self,result,it_atual):
        """
        Metodo modelo material Zener 3D classico, cinematica infinitesimal
        """
        
        # Inicializando variaveis
        F_1=np.zeros((3,3))
        F_0=np.zeros((3,3))
        et_1=np.zeros(6)
        et_0=np.zeros(6)
        tn_1=result.t[it_atual]
        tn_0=result.t[it_atual-1]
        
        sigma_inter=np.zeros(6)
        var_inter=np.zeros(6*self.k_bracos)
        hn_0=result.var_inter[:,it_atual-1]
        Id=np.array([1,1,1,0,0,0])
        cont=0
        
        # delta_t
        dt=np.abs(tn_1-tn_0)
        
        # Retornando tensor F para vetor coluna 
        for i in range(3):
            for j in range(3):
                F_1[j,i]=result.F[cont,it_atual]
                F_0[j,i]=result.F[cont,it_atual-1]
                cont+=1
        
        # Gradiente de deslocamento 
        GradU_1=F_1-np.eye(3)
        GradU_0=F_0-np.eye(3)
        
        # Tensor deformacao infinitesimal
        def_infini_1=0.5*(GradU_1+GradU_1.T)
        def_infini_0=0.5*(GradU_0+GradU_0.T)
        
        # Tensor deformacao infinitesimal em vetor coluna
        for i in range(3):
            et_1[i]=def_infini_1[i,i]
            et_0[i]=def_infini_0[i,i]
        
        et_1[3]=def_infini_1[1,2]; et_1[4]=def_infini_1[0,2]
        et_1[5]=def_infini_1[0,1]; et_0[3]=def_infini_0[1,2]
        et_0[4]=def_infini_0[0,2]; et_0[5]=def_infini_0[0,1]
        
        # Contadores para operar com vetor de variaveis internas
        pa_ini=0
        pa_fim=6
        
        # Traço tensor de deformacao 
        tr_e_1=(et_1[0]+et_1[1]+et_1[2])
        tr_e_0=(et_0[0]+et_0[1]+et_0[2])
        
        # Parcela volumetrica do tensor deformacao
        vol_et_1=(1/3)*tr_e_1*Id
        vol_et_0=(1/3)*tr_e_0*Id
        
        # Parcela desviadora do tensor deformacao 
        dev_et_1=et_1-vol_et_1
        dev_et_0=et_0-vol_et_0
        
        # Parcela desviadora e volumetrica do tensor tensao
        # braço puramente elastico
        S_inf_vol=3*self.K_inf*vol_et_1
        S_inf_dev=2*self.G_inf*dev_et_1
        
        # Parcela viscoelastica desviadora 
        for j in range(self.k_bracos):
            Hhn_0=hn_0[pa_ini:pa_fim]
            A=math.exp(-dt/self.tau[j])*Hhn_0
            B=math.exp(-dt/(2*self.tau[j]))*2*self.G[j]*(dev_et_1-dev_et_0)
            Hhn_1=A+B
            sigma_inter+=Hhn_1
            var_inter[pa_ini:pa_fim]=Hhn_1
            pa_ini+=6
            pa_fim+=6
        
        # Tensao total 
        stress=S_inf_vol+S_inf_dev+sigma_inter
        fail=False
        
        return stress,var_inter,fail,et_1

#-------------------------------------------------------------#     
#--------------------Zener 3D Fracionario---------------------#
#-------------------------------------------------------------# 

class ModeloZener3D_frac_infini(Modelos): 
    """
    Modelo Zener 3D fracionário em cinematica infinitesimal 
    """
    def __init__(self,config):
        """
        Construtor
        """
        self.k_bracos=config['k_bracos']
        self.dimension=6
        self.poisson=config['poisson']
        self.tipo_modelo='fracionario'
#-------------------------------------------------------------#          
    def value_n_time(self,n_time):
        """ 
        Metodo para obter valor final de tempo e inicializar Aj fracionario
        """
        self.n_time=n_time
        self.Aj=np.zeros((self.n_time,self.k_bracos))
#-------------------------------------------------------------#      
    def get_size_var_inter(self):
        """ 
        Metodo para fornecer tamanho vetor variaveis internas
        """
        size=np.zeros(2,dtype=int)
        size[0]=self.dimension*self.k_bracos
        size[1]=self.dimension
        return size

#-------------------------------------------------------------#          
    def set_prop(self,x):
        """ 
        Metodo para atualizar parâmetros do modelo
        """
        self.G=np.zeros(self.k_bracos)
        self.tau=np.zeros(self.k_bracos)
        self.alfa=np.zeros(self.k_bracos)
        self.G_inf=10**x[0]
        cont_1=1
        self.G_0=self.G_inf
        
        for i in range(self.k_bracos):
            self.G[i]=10**x[cont_1+i]
            self.tau[i]=10**x[cont_1+1+i]
            self.alfa[i]=x[cont_1+2+i]
            self.G_0+=self.G[i]
            cont_1+=2
            
        self.K_inf=2*self.G_0*(1+self.poisson)/(3*(1-2*self.poisson))
        
        # Calculo Aj em posicao invertida
        for i in range(self.k_bracos):
            self.Aj[self.n_time-1,i]=1
            cont=self.n_time-2
            j=1
            for k in range(self.n_time-1):
                self.Aj[cont,i]=self.Aj[cont+1,i]*(j-1-self.alfa[i])/j
                cont-=1
                j+=1
        
#-------------------------------------------------------------#         
                
    def mat_solve(self,result,it_atual):
        
        
        """ 
        Metodo modelo material Zener 3D fracionario, cinematica infinitesimal
        """
        F_1=np.zeros((3,3))
        et_1=np.zeros(6)
        tn_1=result.t[it_atual]
        tn=result.t[it_atual-1]
        dt=np.abs(tn_1-tn)
        
        stress_iso_tot=np.zeros(6)
        sigma_inter=np.zeros(6)
        var_inter=np.zeros(6*self.k_bracos)
        cont=0
        
        # Retornando tensor F para vetor coluna
        for i in range(3):
            for j in range (3):
                F_1[j,i]=result.F[cont,it_atual]
                cont+=1
        
        Id=np.array([1,1,1,0,0,0])
        # Gradiente de deslocamento
        GradU_1=F_1-np.eye(3)
        def_infini_1=0.5*(GradU_1+GradU_1.T)
        
        # Tensor deformacao infinitesimal em vetor coluna
        for i in range(3):
            et_1[i]=def_infini_1[i,i]
        
        et_1[3]=def_infini_1[1,2]; et_1[4]=def_infini_1[0,2]
        et_1[5]=def_infini_1[0,1]; 
        
        # Contadores para operar com vetor de variaveis internas
        pa_ini=0
        pa_fim=6
        
        # Traço tensor de deformacao 
        tr_e_1=(et_1[0]+et_1[1]+et_1[2])
                
        # Parcela volumetrica do tensor deformacao
        vol_et_1=(1/3)*tr_e_1*Id
                
        # Parcela desviadora do tensor deformacao 
        dev_et_1=et_1-vol_et_1
                
        # Parcela desviadora e volumetrica do tensor tensao
        # braço puramente elastico
        S_inf_vol=3*self.K_inf*vol_et_1
        S_inf_dev=2*self.G_inf*dev_et_1
        
        
        for j in range(self.k_bracos):
            # Parcela viscoelastica desviadora
            B=1/(self.tau[j]**self.alfa[j])
            A=(dt**(-self.alfa[j])*self.tau[j]**self.alfa[j]+1)*B
            
            # Calculo derivada fracionaria 
            cte_inci=(self.n_time-it_atual-1)
            cte_fim=(self.n_time-1)
            Aj_1=self.Aj[cte_inci:cte_fim,j]
            var_inter_0=result.var_inter[pa_ini:pa_fim,0:it_atual]
            deri_frac_Q=dt**(-self.alfa[j])*np.matmul(var_inter_0,Aj_1)
            # deri_frac_Q=GL_derivative3D(Aj_1,var_inter_0,dt,self.alfa[j])
            
            # Variavel interna iteracao atual 
            var_inter[pa_ini:pa_fim]=(1/A)*(B*2*self.G[j]*dev_et_1-deri_frac_Q)
            
            # Atualizado vetores
            sigma_inter+=var_inter[pa_ini:pa_fim]
            pa_ini+=6
            pa_fim+=6
            
            # Valores elasticos desivadores de cada braco 
            stress_iso_tot+=2*self.G[j]*dev_et_1
            
        stress=S_inf_vol+S_inf_dev+stress_iso_tot-sigma_inter
        fail=False 
        
        return stress,var_inter,fail,et_1
    
#-------------------------------------------------------------#     


    