# -*- coding: utf-8 -*-
"""
Classe abstrata ensaios

@author: mathe
"""
from abc import ABC, abstractmethod
import sys
import numpy as np 

class Ensaios(ABC):
    
    @abstractmethod
    def run(self,material,result):
        """
        Metodo abstrato para rodar ensaio
        """
        pass
  #-------------------------------------------------------------         
  
    @abstractmethod
    def get_n_ensaio(self):
        """
        Metodo abstrato para fornecer n_steps do ensaio 
        """
        pass 
 #-------------------------------------------------------------         
   
    @abstractmethod
    def calc_sse(self,experimento):
        """
        Metodo abstrato para calculo do sse - Otimizacao
        """
        pass 
 #-------------------------------------------------------------         
   
    @abstractmethod
    def get_time(self):
        """
        Metodo abstrato para obter vetor tempo 
        """
        pass
#-------------------------------------------------------------         
   
    def criar_mapeamento_experimentos(self,time_experi):
        """
        Metodo para criar mapeamento entre experimento e numerico
        """
        time_numerico=self.get_time()
        max_time_experi=np.max(time_experi)
        max_time_numerico=np.max(time_numerico)
        
        min_time_experi=np.min(time_experi)
        min_time_numerico=np.min(time_numerico)

        # Verificacoes tamanho vetor numerico e experimental
        if (max_time_experi>max_time_numerico):
            sys.exit('Fatal Error: Tempo final de analise numerica menor que tempo final experimental')
            
        if (min_time_experi<min_time_numerico):
            sys.exit('Fatal Error: Tempo inicial de analise numerica menor que tempo inicial experimental')
        
        
        # Criando mapeamento entre valores experimentais e numericos 
        n_experi=len(time_experi)
        n_numeri=len(time_numerico)
        self.mapeamento_experi=np.zeros(n_experi)
        for i in range(n_experi):
            for j in range(n_numeri-1):
                if (time_experi[i]>=time_numerico[j]) and (time_experi[i]<=time_numerico[j+1]):
                    self.mapeamento_experi[i]=j
                    break 
#-------------------------------------------------------------         
        
        
            
            
        
        

