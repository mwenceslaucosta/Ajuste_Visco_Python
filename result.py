# -*- coding: utf-8 -*-


import numpy as np 

#-------------------------------------------------------------     
class Result:
    
    """
    Classe para inicializar vetores de resultados 
    """
    
    def set_Result(self,material,ensaio):
        size=material.get_size_var_inter()
        size=size.astype(int) #Transformando para inteiro
        size_var_inter=size[0]
        dimension=size[1]
        n_tempo=ensaio.get_n_ensaio()
        self.e=np.zeros((dimension,n_tempo))
        
        # Inicializando F como identidade
        # F=[F11 F21 F31 F12 F22 F32 F13 F23 F33]
        if dimension==6:
            self.F=np.zeros((9,n_tempo))
            self.F[0,:]=1
            self.F[4,:]=1
            self.F[8,:]=1
        else:
            self.F=np.ones(n_tempo)
        
        self.stress=np.zeros((dimension,n_tempo))
        self.t=np.zeros(n_tempo)
        self.var_inter=np.zeros((size_var_inter,n_tempo))
            
 #-------------------------------------------------------------        
        

