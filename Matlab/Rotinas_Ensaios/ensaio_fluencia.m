classdef ensaio_fluencia < Ensaios
    %% --------------------------- Tração ---------------------------------
    %Rotina RUN para ensaio de fluência.
    
    properties
        
        sigma_f             %Tensão de fluência
        t_f                 %Tempo de fluência
        deltat              %intervalo de tempo
        t                   %vetor tempo discretizado
        n                   %Número de intervalos de análise
        e                   %Historico de deformação
        N_iter              %Numero mmaximo de iterações
        eps                 %Variação para dif. finitas
        tol_conv            %Erro para criterio de convergencia newton
        t_carregamento
        
    end
    
    methods
        %% Construtor ensaio de Fluência
        function obj = ensaio_fluencia()
        end
        
        function set_ensaio_fluencia(obj,fluencia)
            obj.sigma_f =fluencia.stress;
            obj.t_f=fluencia.tf;
            obj.deltat=fluencia.deltat;    %Intervalo de tempo [s]
            obj.t=(0:obj.deltat:obj.t_f);  %tempo discretizado
            obj.t_carregamento=fluencia.t_carregamento; %Controla se será utilizada rampa de carregamento ou aplicacao instantanea
            obj.n=size(obj.t,2);
            obj.N_iter=100;
            obj.eps=1E-6;
            obj.tol_conv=1E-6;
        end
        %% Método que fornece tamanho vetor t para inicializar variaveis.
        function n=get_n_ensaio(obj)
            n=obj.n;
        end
        %% Método para fornecer vetor tempo numérico
        function time=get_time(obj)
            time=obj.t;
        end
        
        %% Método construção rampa de tração para fluência
        
        function [rampa_stress_e_defor]=rampa_stress_hist_def(obj,material)
            size=material.get_size_var_inter();
            dimension=size(2);
            rampa_stress_e_defor.dimension=dimension;
            
            %Rampa stress
            rampa_stress_e_defor.hist_stress=zeros(dimension,obj.n);
            %obj.t_carregamento(1)=1 -> Rampa linear de carregamento em Fluencia
            %obj.t_carregamento(1)=2 -> Carregamento instantaneo em Fluencia
            %obj.t_carregamento(1)=3 -> Rampa linear de carregamento e
            %descarregamento em Fluencia seguido por recuperacao.
            
            %obj.t_carregamento(2)-> Tempo de carregamento (s) ate tensao de
            %fluencia.
            switch obj.t_carregamento(1)
                
                case 1 %Rampa de carregamento linear - Somente Fluencia
                    a=obj.sigma_f/obj.t_carregamento(2);
                    for i=1:obj.n
                        if obj.t(i)<=obj.t_carregamento(2)
                            rampa_stress_e_defor.hist_stress(1,i)=a*obj.t(i);
                        else
                            rampa_stress_e_defor.hist_stress(1,i)=obj.sigma_f;
                        end
                    end
                   
                case 2  %Carregamento instantaneo - Somente Fluencia 
                    rampa_stress_e_defor.hist_stress(1,2:obj.n)=obj.sigma_f;
                    
                    %% Fluencia e recuperacao - Rampa Linear 
                    %Rampa de carregamento, seguido por fluencia em tensao
                    %sigma_f, seguido por descarregamento linear e recuperacao em
                    %tensao nula durante tempo de recuperacao obj.t_carregamento(3).
                   
                    %obj.t_carregamento(2)-> Tempo de carregamento(s) e
                    %descarregamento.
                    %obj.t_carregamento(3)-> Tempo de fluencia e recuperacao
                case 3
                    t_1=obj.t_carregamento(2);
                    t_2=obj.t_carregamento(2)+obj.t_carregamento(3);
                    t_3=2*obj.t_carregamento(2)+obj.t_carregamento(3);
                    a_1=obj.sigma_f/t_1;
                    a_2=(obj.sigma_f*(1-1/(1-t_2/t_3)))/t_2;
                    b_2=obj.sigma_f/(1-t_2/t_3);
                    for i=1:obj.n
                        if obj.t(i)<=t_1
                            rampa_stress_e_defor.hist_stress(1,i)=a_1*obj.t(i);
                        elseif obj.t(i)>=t_1 && obj.t(i)<=t_2
                            rampa_stress_e_defor.hist_stress(1,i)=obj.sigma_f;
                        elseif obj.t(i)>=t_2 && obj.t(i)<=t_3
                            rampa_stress_e_defor.hist_stress(1,i)=a_2*obj.t(i)+b_2;
                        else 
                            rampa_stress_e_defor.hist_stress(1,i)=0.07;
                        end
                    end
            end
        end
        
        %% ---------------------RUN Fluencia --------------------------%%
        
        function [result,fail]=run(obj,material,result)
            
            [rampa_stress_e_defor]=rampa_stress_hist_def(obj,material);
            result.t=obj.t;
            
            %% Condicionamento para chamar equilibrio para cinematica
            
            for i=2:obj.n
                %Chamada newton para solução do equilibrio mecânico
                [stress,var_inter,fail,result,strain] = equilibrio_fluencia(obj,...
                    material,result,i,rampa_stress_e_defor) ;
                %Fim newton passo i
                
                if fail %Se verdadeiro, falha e retorna
                    return
                end
                
                result.stress(:,i)=stress;
                result.e(:,i)=strain;
                result.var_inter(:,i)=var_inter;
            end
            
            %% Newton para equilibrio mecanico cinematica linearizada
            function [stress,var_inter,fail,result,strain]=equilibrio_fluencia(obj,...
                material,result,it_atual,rampa_stress_e_defor)
                hist_stress=rampa_stress_e_defor.hist_stress;
                
                %Atualizando com anteriores. 
                result.F(:,it_atual)=result.F(:,it_atual-1);
                
                for j=1:obj.N_iter
                    
                    %Verificação residuo
                    [stress,var_inter,fail,strain]=material.mat_solve(result,it_atual);
                    
                    dimension=rampa_stress_e_defor.dimension;
                    if dimension ==6
                        n_var_residuo=2;
                        n_componentes_perturbadas_F=[1 5];
                    else
                        n_var_residuo=1;
                        n_componentes_perturbadas_F=1;
                    end
                    R=stress(1:n_var_residuo,1)-hist_stress(1:n_var_residuo,it_atual);
                    if fail
                        return
                    end
                    
                    norm_R=norm(R);
                    if norm_R < obj.tol_conv
                        return
                    end
                    %% modulo tangente numerico
                    F_0 = result.F(:,it_atual);
                    Dif=zeros(n_var_residuo);
                    %% Diferenças finitas centrais
                    for c = 1:n_var_residuo
                        n_pertubacao=n_componentes_perturbadas_F(c);
                        result.F(n_pertubacao , it_atual ) = F_0(n_pertubacao) + obj.eps;
                        if n_var_residuo==2
                            result.F(9 , it_atual ) = result.F(5 , it_atual );
                        end
                        [stress_front,~,fail]=material.mat_solve(result,it_atual);
                        if fail
                            return
                        end
                        
                        result.F(n_pertubacao,it_atual) = F_0(n_pertubacao) - obj.eps ;
                        if n_var_residuo==2
                            result.F(9 , it_atual ) = result.F(5 , it_atual );
                        end
                        [stress_back,~,fail]=material.mat_solve(result,it_atual);
                        if fail
                            return
                        end
                        Dif(:,c) = (stress_front(1:n_var_residuo) - stress_back(1:n_var_residuo)) / (2*obj.eps) ;
                        result.F(n_pertubacao,it_atual) = F_0(n_pertubacao); %Retornando vetor deformação para condição do passo i anterior a dif. finitas
                         if n_var_residuo==2
                            result.F(9 , it_atual ) = result.F(5 , it_atual );
                         end
                    end
                    
                    if det(Dif)==0
                        fail=true;
                        disp('singular, retornou');
                        return
                    end
                    
                    delta_F=-linsolve(Dif,R);
                    for k=1:n_var_residuo
                    n_pertubacao=n_componentes_perturbadas_F(k);
                    result.F(n_pertubacao,it_atual)=result.F(n_pertubacao,it_atual)+delta_F(k);
                    end 
                    if n_var_residuo==2
                        result.F(9 , it_atual ) = result.F(5 , it_atual );
                    end
                    
                end %Fim iterações Newton
                fail=true; %Se chegar aqui, newton não convergiu para 100 iter
            end %Fim equilibrio mecanico fluencia cinematica linearizada
            
            
        end %Fim run fluencia
        
        %% Método para calculo do SSE ensaio de fluencia
        
        function sse=calc_sse(obj,experimento,result)
            sse=0;
            for i=1:length(experimento.t)
                i_numerical=obj.mapeamento_experi(i); %Indice numerico correspondente
                factor=(experimento.t(i)-result.t(i_numerical))/(result.t(i_numerical+1)-result.t(i_numerical));
                numerical_value=result.e(1,i_numerical)+factor*(result.e(1,i_numerical+1)-result.e(1,i_numerical));
%               sse=(experimento.e(i,1)-numerical_value)^2*experimento.e(i,2)+sse;
                sse=(experimento.e(i,1)-numerical_value)^2+sse;
            end
%               sse=sse/sum(experimento.e(i,2));

          %% Nao precisa de valores interpolados - utiliza vetor tempo corrigido
%             for i=1:length(experimento.t)
%                 i_numerical=obj.mapeamento_experi(i); %Indice numerico correspondente
%                 numerical_value=result.e(1,i_numerical);
%                 sse=(experimento.e(i,1)-numerical_value)^2+sse;
%             end
            
%             for i=1:length(experimento.t)
%                 i_numerical=obj.mapeamento_experi(i); %Indice numerico correspondente
%                 factor=(experimento.t(i)-result.t(i_numerical))/(result.t(i_numerical+1)-result.t(i_numerical));
%                 numerical_value=result.e(1,i_numerical)+factor*(result.e(1,i_numerical+1)-result.e(1,i_numerical));
%                 sse=(norm(experimento.e(i)-numerical_value)/experimento.e(i))^2+sse;
%             end
%                sse=sqrt(sse);
        end %Fim método calculo SSE
        
    end
end