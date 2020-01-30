# -*- coding: utf-8 -*-
"""
Modulo contendo Ensaios numericos: Tracao, fluencia-recuperacao, cisalhamento

@author: mathe
"""
from ensaios import Ensaios 
import numpy as np
from scipy import linalg 

#-------------------------------------------------------------# 
#-------------------------------------------------------------# 
#-------------------------------------------------------------#     
#--------------------Ensaio de Tracao-------------------------#
#-------------------------------------------------------------#
#-------------------------------------------------------------# 
#-------------------------------------------------------------# 
class Ensaio_tracao(Ensaios):
    """
    Classe ensaio de tracao
    """
    
    def set_ensaio(self,tensile):
        """
        Metodo para setar ensaio de tracao
        """
        self.def_rate=tensile['def_rate']
        self.def_max=tensile['def_max']
        self.tens_tf=self.def_max/self.def_rate
        self.step_tempo=tensile['step_tempo']
        self.t=np.linspace(0,self.tens_tf,self.step_tempo)
        self.n=self.step_tempo
        self.N_iter=100
        self.eps=1E-6
        self.tol_conv=1E-6
        self.tens_deltat=self.tens_tf/self.step_tempo
    
#-------------------------------------------------------------#         

    def get_n_ensaio(self):
        """
        Metodo para fornecer tamanho do vetor tempo
        """
        n=self.n
        return n
#-------------------------------------------------------------#  
    def get_time(self):
        """
        Metodo para fornecer vetor tempo numerico
        """
        time=self.t
        return time 
         
#-------------------------------------------------------------#      
    def historico_def(self,material):
        """
        Metodo para fornecer vetor historico de deformacao numerico
        """
        strain_hist={}
        size=material.get_size_var_inter()
        size=size.astype(int)
        dimension=size[1]
        strain_hist['dimension']=dimension
                
        # Historico deformacao 
        strain_hist['hist_strain']=np.zeros((dimension,self.n))
        strain_hist['hist_strain'][0,:]=self.def_rate*self.t
        
        
        # Recuperando F pela deformacao de engenharia 
        strain_hist['F_11']=strain_hist['hist_strain'][0,:]+1
        return strain_hist
#-------------------------------------------------------------# 
    
    def run(self,material,result):
        
        """
        Metodo run ensaio de tracao
        """
        strain_hist=self.historico_def(material)
        result.e=strain_hist['hist_strain']
        result.F[0,:]=strain_hist['F_11']
        result.t=self.t
     
        # Loop RUN
        for i in range(1,self.n):
            # Equilibrio caso 3D
            if (strain_hist['dimension']==6):
                # Chamada Newton para equilibrio mecanico do ensaio
                stress,var_inter,fail,result=self.equilibrio_tracao(material,result,i,strain_hist)
            else: 
                stress,var_inter,fail,result=material.mat_solve(result,i)
            
            if fail: 
                return result,fail
                         
            result.stress[:,i]=stress
            result.var_inter[:,i]=var_inter
                              
        return result,fail 

#-------------------------------------------------------------#     
    def equilibrio_tracao(self,material,result,it_atual,strain_hist):
        """
        Newton para equilibrio mecanico do ensaio de tracao 
        """
        # Atualizando com componentes do passo anterior.
        # Exceto componente 0. 
        result.F[1:9,it_atual]=result.F[1:9,it_atual-1]
        
        for i in range(self.N_iter):
            # Chamada modelo material 
            stress,var_inter,fail,et_1=material.mat_solve(result,it_atual)
            result.e[1:6,it_atual]=et_1[1:6]
            
            if fail:
                return stress,var_inter,fail,result
            
            # Verificacao residuo 
            R=stress[1]
            # Norma do residuo 
            norm_R=linalg.norm(R)
            if norm_R<self.tol_conv:
                return stress,var_inter,fail,result
            
            # Modulo Tangente numerico 
            F_0=np.zeros(9)
            F_0[:]=result.F[:,it_atual] # Atencao em objetos mutaveis
            # Diferencas finitas centrais 
            # Lembre-se que componente [0,0] é dada. 
            
            result.F[4,it_atual]=F_0[4]+self.eps 
            result.F[8,it_atual]=result.F[4,it_atual]
            stress_front,var_inter,fail,et_1=material.mat_solve(result,it_atual)
            
            if fail:
                return stress,var_inter,fail,result
                        
            result.F[4,it_atual]=F_0[4]-self.eps 
            result.F[8,it_atual]=result.F[4,it_atual]
            stress_back,var_inter,fail,et_1=material.mat_solve(result,it_atual)
            
            if fail:                 
                return stress, var_inter_,fail,result
            
            # Obtendo solucao do sistema linear e obtendo incremento para F.
            Dif=(stress_front-stress_back)/(2*self.eps)
            result.F[4,it_atual]=F_0[4]
            
            delta_F=-R/Dif[1]
            result.F[4,it_atual]=result.F[4,it_atual]+delta_F
            result.F[8,it_atual]=result.F[4,it_atual]
            # Fim iteracoes newton       
        
        print('Newton ensaio de traçao nao convergiu')
        fail=True
        
        return stress, var_inter,fail,result
    
#-------------------------------------------------------------#  
    
    def  calc_sse(self,experimento,result):
        """
        Metodo para calculo do SSE ensaio de traçao
        """
        n=len(experimento.t)
        sse=0
        for i in range(n):
            i_numerical=self.mapeamento_experi[i]
            A=experimento.t[i]-result.t[i_numerical]
            B=result.t[i_numerical+1]-result.t[i_numerical]
            factor=A/B
            C=factor*(result.stress[0,i_numerical+1]-result.stress[0,i_numerical])
            numerical_value=result.stress[0,i_numerical]+C
            sse+=(experimento.stress[i]-numerical_value)**2
        
        return sse
        
#-------------------------------------------------------------#        
#-------------------------------------------------------------#       
#-------------------------------------------------------------#        
#------------Ensaio de fluencia-recuperacao-------------------#           
#-------------------------------------------------------------#    
#-------------------------------------------------------------#       
#-------------------------------------------------------------#
        
class Ensaio_fluencia_recuperacao(Ensaios):
    """
    Classe ensaio de fluencia e fluencia-recuperação
    """
        
#-------------------------------------------------------------# 
        
    def set_ensaio(self,fluencia):
        """
        Metodo para setar ensaio de fluencia
        """
        self.sigma_f=fluencia['stress'] #Tensao de fluencia
        self.t_f=fluencia['tf']         #tempo final de fluencia
        self.deltat=fluencia['delta_t']  #Intervalo de tempo numerico 
        self.t_carregamento=fluencia['t_carregamento'] #Vetor com dados de
                                                           #carregamento
        self.t=np.arange(0,self.t_f+self.deltat,self.deltat) 
        self.n=len(self.t)
        self.N_iter=100
        self.eps=1E-6
        self.tol_conv=1E-6
        if self.t_carregamento[0]==3:
            self.sigma_recovery=fluencia['sigma_recovery']
#-------------------------------------------------------------#  
        
    def get_n_ensaio(self):
        """
        Metodo para fornecer tamanho do vetor t 
        """            
        n=self.n
        return n
#-------------------------------------------------------------#  
            
    def get_time(self):
        """ 
        Metodo para fornecer vetor tempo numerico
        """
        time=self.t
        return time 
#-------------------------------------------------------------#  
        
    def rampa_stress_hist_def(self,material):
        """
        Metodo para criar rampa de fluencia ou fluencia-recuperacao 
        """
        size=material.get_size_var_inter()
        dimension=size[1]
        rampa_stress_e_defor={}
        rampa_stress_e_defor['dimension']=dimension
            
        #Criando rampa de carregamento e descarregamento 
        #self.t_carregamento[0]: Seleciona tipo de ensaio de fluencia
        #self.t_carregamento[0]=1 -> Rampa linear de carregamento em Fluencia
        #self.t_carregamento[0]=2 -> Carregamento instantaneo em Fluencia
        #self.t_carregamento[0]=3 -> Rampa linear de carregamento e
        #descarregamento em Fluencia seguido por recuperacao.
            
        #self.t_carregamento(2):Tempo de carregamento (s) ate tensao de
        #fluencia
            
        rampa_stress_e_defor['hist_stress']=np.zeros((dimension,self.n))
          
        #Caso 1 - Rampa de carregamento linear - Somente fluencia
        if (self.t_carregamento[0] == 1):               
            a=self.sigma_f/self.t_carregamento[1]
            for i in range(self.n):
                if (self.t[i]<=self.t_carregamento[1]):
                    rampa_stress_e_defor['hist_stress'][0,i]=a*self.t[i]
                else:
                    rampa_stress_e_defor['hist_stress'][0,i]=self.sigma_f
            
        #Caso 2 - Carregamento instantaneo - Somente fluencia
        if (self.t_carregamento[0] == 2):                
            rampa_stress_e_defor['hist_stress'][0,1:(self.n+1)]=self.sigma_f
            
        #Caso 3 - Fluencia e recuperacao - Rampa Linear no carregamento e 
        #descarregamento.
            
        #Rampa de carregamento linear durante tempo self.t_carregamento[1], 
        #seguido por fluencia em tensao sigma_f durante tempo
        # self.t_carregamento[2], seguido por descarregamento linear 
        # durante tempo self.t_carregamento[1] e recuperacao em tensao nula 
        #durante tempo de recuperacao self.t_carregamento[2].  
        if (self.t_carregamento[0] == 3): 
            t_1=self.t_carregamento[1] 
            #t_1: Tempo de carregamento 
            t_2=self.t_carregamento[1]+self.t_carregamento[2] 
            #t_2: Tempo de carregamento + Fluencia
            t_3=2*self.t_carregamento[1]+self.t_carregamento[2] 
            #t_3: Tempo de inicio da recuperacao
            
            a_1=self.sigma_f/t_1
            a_2=(self.sigma_f*(1-1/(1-t_2/t_3)))/t_2
            b_2=self.sigma_f/(1-t_2/t_3);
                
            for i in range(self.n):
                if (self.t[i]<=t_1):
                    rampa_stress_e_defor['hist_stress'][0,i]=a_1*self.t[i]
                elif (self.t[i]>=t_1 and self.t[i]<=t_2):
                    rampa_stress_e_defor['hist_stress'][0,i]=self.sigma_f
                elif (self.t[i]>=t_2 and self.t[i]<=t_3):
                    rampa_stress_e_defor['hist_stress'][0,i]=a_2*self.t[i]+b_2
                else:
                    rampa_stress_e_defor['hist_stress'][0,i]=self.sigma_recovery
            
        return rampa_stress_e_defor
 
 #-------------------------------------------------------------#               

    def run(self,material,result):
        
        """
        Metodo run ensaio de fluencia e fluencia-recuperacao
        """
        rampa_stress_e_defor=self.rampa_stress_hist_def(material)
        result.t=self.t
         
        # Loop RUN
        for i in range(1,self.n):
            # Equilibrio mecanico                
            # Chamada Newton para equilibrio mecanico do ensaio
            stress,var_inter,fail,result,strain=self.equilibrio_fluencia(material,result,i,rampa_stress_e_defor)                
                
            if fail: 
                return result,fail
                             
            result.stress[:,i]=stress
            result.e[:,i]=strain
            result.var_inter[:,i]=var_inter
                                  
        return result,fail    

#-------------------------------------------------------------#                     
        
    def equilibrio_fluencia(self,material,result,it_atual,rampa_stress_e_defor):
        """
        Metodo com Newton para equilibrio mecanico em fluencia
        """
        hist_stress=rampa_stress_e_defor['hist_stress']
        dimension=rampa_stress_e_defor['dimension']
            
        #Definindo quais componentes de F sao perturbadas. 
        if (dimension==6):
            n_var_residuo=2
            n_componentes_perturbadas_F=np.array([0,4])
        else:
            n_var_residuo=1
            n_componentes_perturbadas_F=0
                    
        #Atualizando com iteracoes anteriores
        result.F[:,it_atual]=result.F[:,it_atual-1]
            
        for j in range(self.N_iter):
                
            # Chamada modelo material 
            stress,var_inter,fail,strain=material.mat_solve(result,it_atual)
            if fail:                    
                return stress,var_inter,fail,result,strain
                
            #Verificacao residuo 
            R=stress[0:n_var_residuo]-hist_stress[0:n_var_residuo,it_atual]
            norm_R=norm_R=linalg.norm(R)
            if (norm_R < self.tol_conv):
                return stress,var_inter,fail,result,strain
                
            #Modulo tangente numerico 
            F_0=np.zeros(len(result.F[:,it_atual]))
            F_0[:]=result.F[:,it_atual]
            Dif=np.zeros((n_var_residuo,n_var_residuo))
                
            #Diferenças finitas centrais 
            for k in range(n_var_residuo):
                n_pertubacao=n_componentes_perturbadas_F[k]
                #Perturbaçao a frente 
                result.F[n_pertubacao,it_atual]=F_0[n_pertubacao]+self.eps
                if (n_var_residuo==2):
                    result.F[8,it_atual]=result.F[4,it_atual]
                    
                stress_front,var_inter,fail,et_1=material.mat_solve(result,it_atual)
                if fail:
                    return stress,var_inter,fail,result,strain
                    
                #Perturbaçao atras
                result.F[n_pertubacao,it_atual]=F_0[n_pertubacao]-self.eps
                if (n_var_residuo==2):
                    result.F[8,it_atual]=result.F[4,it_atual]
                    
                stress_back,var_inter,fail,et_1=material.mat_solve(result,it_atual)
                if fail:
                    return stress,var_inter,fail,result,strain
                    
                Dif[:,k]=(stress_front[0:n_var_residuo]-stress_back[0:n_var_residuo])/(2*self.eps)
                result.F[n_pertubacao,it_atual]=F_0[n_pertubacao]
                if (n_var_residuo==2):
                    result.F[8,it_atual]=result.F[4,it_atual]
                
            #Obtendo incremento para F
            delta_F=-linalg.solve(Dif,R)
            for k in range(n_var_residuo):
                n_pertubacao=n_componentes_perturbadas_F[k]
                result.F[n_pertubacao,it_atual]=result.F[n_pertubacao,it_atual]+delta_F[k]
                
            if (n_var_residuo==2):
                result.F[8,it_atual]=result.F[4,it_atual]
            
        print('Newton ensaio de fluencia nao convergiu')            
        fail=True  #Fim iteracoes Newton - Se chegar aqui falhou.       
        return stress,var_inter,fail,result,strain 

#-------------------------------------------------------------#                     
    
    def  calc_sse(self,experimento,result):
        """
        Metodo para calculo do SSE ensaio de fluencia
        """
        n=len(experimento.t)
        sse=0
        for i in range(n):
            i_numerical=self.mapeamento_experi[i].astype(int)
            A=experimento.t[i]-result.t[i_numerical]
            B=result.t[i_numerical+1]-result.t[i_numerical]
            factor=A/B
            C=factor*(result.e[0,i_numerical+1]-result.e[0,i_numerical])
            numerical_value=result.e[0,i_numerical]+C
            sse+=(experimento.e[i]-numerical_value)**2
        
        return sse
        
            

