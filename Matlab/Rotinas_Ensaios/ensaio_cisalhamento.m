classdef ensaio_cisalhamento < Ensaios
    %% ------------------------- Cisalhamento -------------------------
    %Rotina RUN para ensaio de tra��o.
    %% ------------- Informa��es variaveis ensaio de tra��o --------------- %%
    %tensile.rate->Taxa de deforma��o p/dire��o 1, 2 e 3
    %tensile.tf->   Tempo final em segundos
    %tensile.deltat-> Incremento de tempo
    %recebe deforma��o [epslon12]
    
    %% --------------- Objetos particulares desta classe -------------- %%
    properties
        t                %vetor tempo discretizado
        e                %Hist�rico de deforma��o
        n                %N�mero de intervalos de an�lise
        def_rate
        def_max
        step_tempo
        shear_tf
        shear_deltat
        
    end
    
    methods
        %% ---------------- Construtor Ensaio de tra��o---------------- %%
        function obj= ensaio_cisalhamento()
        end 
        
        function  set_ensaio_cisalhamento(obj,shear)
            obj.def_rate=shear.def_rate;                   %Taxa def. dire��o 1
            obj.def_max=shear.def_max;                     %Maxima deforma��o ensaio
            obj.shear_tf=obj.def_max/obj.def_rate;         %Tempo final ensaio
            obj.step_tempo=shear.step_tempo;               %Intervalo de tempo [s]
            obj.t=linspace(0,obj.shear_tf,obj.step_tempo); %Vetor tempo discretizado
            obj.n=size(obj.t,2);
            obj.shear_deltat=obj.shear_tf/obj.step_tempo;
        end
        
        %% M�todo que fornece tamanho vetor t para inicializar variaveis.
        function n=get_n_ensaio(obj)
            n=obj.n;
        end
        %% M�todo para fornecer vetor tempo num�rico
        function time=get_time(obj)
            time=obj.t;
        end
        
        %% Historico de deforma��o 
        
        function  strain_hist=historico_def(obj,material)
                  size=material.get_size_var_inter();
                  dimension=size(2);
                  strain_hist.dimension=dimension;
                  
                  %Historico
                  strain_hist.hist_strain=zeros(dimension,obj.n);
                  %F=[F11 F21 F31 F12 F22 F32 F13 F23 F33];
                  strain_hist.hist_strain(6,:)=obj.def_rate*obj.t;
                  strain_hist.F_12=strain_hist.hist_strain(6,:)*2;
        end 
        %% ------------------ Run Ensaio Cisalhamento------------------- %%
        
        function [result,fail]=run(obj,material,result)
            
            strain_hist=historico_def(obj,material);
            result.e=strain_hist.hist_strain;
            result.t=obj.t;
            result.F(4,:)=strain_hist.F_12; %Alocando posicao 12 de F
            
            for i=2:obj.n
                [stress,var_inter,fail]=material.mat_solve(result,i);   
                if fail %Se verdadeiro, falha e retorna
                    return
                end
                
                result.stress(:,i)=stress;
                result.var_inter(:,i)=var_inter;
                 
            end
        end %Fim Run ensaio de tra��o
        
        %% M�todo para calculo do SSE ensaio de tra��o
        
        function sse=calc_sse(obj,experimento,result)
            sse=0;
            for i=1:length(experimento.t)
                i_numerical=obj.mapeamento_experi(i); 
                factor=(experimento.t(i)-result.t(i_numerical))/(result.t(i_numerical+1)-result.t(i_numerical));
                numerical_value=result.stress(6,i_numerical)+factor*(result.stress(6,i_numerical+1)-result.stress(6,i_numerical));
                sse=(experimento.stress(i)-numerical_value)^2+sse;
            end
        end %Fim m�todo calculo SSE
    end %Fim metodos
end %Fim classe ensaios de tra��o
