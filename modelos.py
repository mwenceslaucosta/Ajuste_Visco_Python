# -*- coding: utf-8 -*-
"""
Classe abstrata modelos

@author: mathe
"""

from abc import ABC, abstractmethod

class Modelos(ABC):
    """
    Classe abstrata Modelos 
   
    """
    @abstractmethod
    def mat_solve(self,result,i):
         """
         Metodo abstrato que fornece resultado do modelo material 
         """
         pass
         #print(modelos.Modelos.mat_solve.__doc__)
    @abstractmethod
    def get_size_var_inter(self):
        """
        Metodo abstrato que fornece tamanho do vetor de var_internas    do modelo material 

        """
        pass 
   
    @abstractmethod
    def set_prop(self,x):
        """
        Metodo abstrato que atribui valor as variaiveis do modelo 
        """
        pass
        
     

 
    



        
