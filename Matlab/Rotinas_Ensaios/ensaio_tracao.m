classdef ensaio_tracao < Ensaios
    %% --------------------------- Tração ---------------------------------
    %Rotina RUN para ensaio de tração.
    %% ------------- Informações variaveis ensaio de tração --------------- %%
    %tensile.rate->Taxa de deformação p/direção 1, 2 e 3
    %tensile.tf->   Tempo final em segundos
    %tensile.deltat-> Incremento de tempo
    
    %% --------------- Objetos particulares desta classe -------------- %%
    properties
        tens_tf          %Tempo final
        tens_deltat      %intervalo de tempo
        t                %vetor tempo discretizado
        e                %Histórico de deformação
        n                %Número de intervalos de análise
        N_iter           %Numero mmaximo de iterações
        eps              %Variação para dif. finitas
        tol_conv         %Erro para criterio de convergencia newton
        def_rate
        def_max
        step_tempo
        
    end
    
    methods
        %% ---------------- Construtor Ensaio de tração---------------- %%
        function obj=ensaio_tracao()
        end 
        
        function  set_ensaio_tracao(obj,tensile)
            obj.def_rate=tensile.def_rate;                 %Taxa def. direção 1
            obj.def_max=tensile.def_max;                   %Maxima deformação ensaio
            obj.tens_tf=obj.def_max/obj.def_rate;          %Tempo final ensaio
            obj.step_tempo=tensile.step_tempo;             %Intervalo de tempo [s]
            obj.t=linspace(0,obj.tens_tf,obj.step_tempo);  %Vetor tempo discretizado
            obj.n=size(obj.t,2);
            obj.N_iter=100;
            obj.eps=1E-6;
            obj.tol_conv=1E-6;
            obj.tens_deltat=obj.tens_tf/obj.step_tempo;
        end
        
        %% Método que fornece tamanho vetor t para inicializar variaveis.
        function n=get_n_ensaio(obj)
            n=obj.n;
        end
        %% Método para fornecer vetor tempo numérico
        function time=get_time(obj)
            time=obj.t;
        end
        
        %% Historico de deformação 
        
        function  strain_hist=historico_def(obj,material)
                  size=material.get_size_var_inter();
                  dimension=size(2);
                  strain_hist.dimension=dimension;
                  
                  %Historico
                  strain_hist.hist_strain=zeros(dimension,obj.n);
                  strain_hist.hist_strain(1,:)=obj.def_rate*obj.t;
                  %Recuperando pela deformacao verdadeira. %DISCUTIR - JAN
                  %F=[F11 F21 F31 F12 F22 F32 F13 F23 F33];
%                 strain_hist.F_11=exp(strain_hist.hist_strain(1,:)); 
                  %Recuperando pela deformacao de engenharia. 
                  strain_hist.F_11=strain_hist.hist_strain(1,:)+1;
        end 
        %% ------------------ Run Ensaio Tração----------------------- %%
        
        function [result,fail]=run(obj,material,result)
            
            strain_hist=historico_def(obj,material);
            result.e=strain_hist.hist_strain;
            result.F(1,:)=strain_hist.F_11;
            result.t=obj.t;
            
            for i=2:obj.n                
                 if strain_hist.dimension==6 %Entra somente se for 3D
                %Chamada newton para solução do equilibrio mecânico
                [stress,var_inter,fail,result,strain] = equilibrio_tracao(obj,material,result,i,strain_hist) ;
                 else
                [stress,var_inter,fail,strain]=material.mat_solve(result,i);
                 end
                %Fim newton passo i                 
                if fail %Se verdadeiro, falha e retorna
                    return
                end                
                result.stress(:,i)=stress;
                result.e(:,i)=strain;
                result.var_inter(:,i)=var_inter;                 
            end
            
            %% Newton para equilibrio mecanico 
            function [stress,var_inter,fail,result,strain]=equilibrio_tracao(obj,material,result,it_atual,strain_hist)

                %Atualizando com anteriores. Exceto para componente 1
                result.F(2:9,it_atual)=result.F(2:9,it_atual-1);
                for j=1:obj.N_iter
                    
                    %Chamada modelo constitutivo 
                    [stress,var_inter,fail,strain]=material.mat_solve(result,it_atual);
                    
                    if fail
                        return
                    end
                   %Verificação residuo (3D) - Simplificado para modelo
                   %isotropico em tracao
                    R=stress(2); 
                    norm_R=norm(R);
                    if norm_R < obj.tol_conv
                        return
                    end
                    %% modulo tangente numerico
                    F_0=result.F(:,it_atual);
                    %% Diferenças finitas centrais
                    %Obs.: Componente (1,1) é dada                   
                        result.F(5 , it_atual ) = F_0(5) + obj.eps;
                        result.F(9 , it_atual ) = result.F(5 , it_atual );
                        [stress_front,~,fail]=material.mat_solve(result,it_atual);
                        if fail
                            return
                        end
                        result.F(5,it_atual) = F_0(5) - obj.eps ;
                        result.F(9 , it_atual ) = result.F(5 , it_atual );
                        [stress_back,~,fail]=material.mat_solve(result,it_atual);
                        if fail
                            return
                        end
                        Dif = (stress_front - stress_back) / (2*obj.eps) ; 
                        result.F(5,it_atual) = F_0(5); %Retornando vetor deformação para condição do passo i anterior a dif. finitas
                    
                    delta_F=-R/Dif(2); 
                    result.F(5,it_atual)=result.F(5,it_atual)+delta_F; 
                    result.F(9 , it_atual )= result.F(5 , it_atual );
                end %Fim iterações Newton
                fail=true; %Se chegar aqui, newton não convergiu para 100 iter
            end %Fim Newton para equilibrio mecanico
        end %Fim Run ensaio de tração
        
        %% Método para calculo do SSE ensaio de tração
        
        function sse=calc_sse(obj,experimento,result)
            sse=0;
            for i=1:length(experimento.t)
                i_numerical=obj.mapeamento_experi(i); 
                factor=(experimento.t(i)-result.t(i_numerical))/(result.t(i_numerical+1)-result.t(i_numerical));
                numerical_value=result.stress(1,i_numerical)+factor*(result.stress(1,i_numerical+1)-result.stress(1,i_numerical));
                sse=(experimento.stress(i)-numerical_value)^2+sse;
            end
        end %Fim método calculo SSE
    end %Fim metodos
end %Fim classe ensaios de tração

