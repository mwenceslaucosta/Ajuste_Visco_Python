classdef Result < handle
    %% Classe para inicializar vetores de resultado
    
    %% ----------------------Informações Preparação mémoria -----------------%%
    %Seção para inicializar deformação, tensão e variáveis internas.
    %size_var_inter -> Número de variaveis internas do modelo em análise
    %n_tempo -> Número de steps de tempo
    %result.e -> Deformação  ; result.stress -> Tensão
    %result.t -> tempo       ; result.var_inter -> Variaveis internas
    %%
    properties
        t
        e
        stress
        var_inter
        F
    end
    
    methods
        
        function set_Result(obj,material,ensaio)
            size=material.get_size_var_inter();
            size_var_inter=size(1);
            dimension=size(2);
            n_tempo=ensaio.get_n_ensaio();
            obj.e=zeros(dimension,n_tempo);
            %Iniciallizando F como identidade
            %F=[F11 F21 F31 F12 F22 F32 F13 F23 F33];
            if dimension==6
            obj.F=zeros(9,n_tempo);
            obj.F(1,:)=1; 
            obj.F(5,:)=1; obj.F(9,:)=1; 
            else 
            obj.F=ones(1,n_tempo);
            end 
            obj.stress=zeros(dimension,n_tempo);
            obj.t=zeros(1,n_tempo);
            obj.var_inter=zeros(size_var_inter,n_tempo);
            
        end
        
    end
    
end
