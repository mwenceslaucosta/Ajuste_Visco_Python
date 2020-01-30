classdef zener3D_fracionario_funcao_de_G_e_K_em_tensao<modelos
    %% Modelo Viscoel�stico 3D - cinem�tica linearizada
    %% Com parcela viscoelastica somente desviadora
    %% ----------Informa��es modelo zener 3D fracionario --------------- %%
    %k_bracos-> Numero de bra�os. Valor n�o atualizado durante otimiza��o
    %x-> Vet. de projeto contendo par�metros materiais.
    % x(1)=G_inf; x(2)=G_1; x(3)=tau1; x(4)=alfa_1 ;x(5)=G_2;  x(6)=tau2  x(7)=alfa2 ...... 
    
    %fornece deforma��es [e11 e22 e33 e23 e13 e12]
    
    properties
        G
        G_inf
        K_inf
        tau
        alfa
        k_bracos
        dimension
        poisson
        flag_cinematica
        n_time
        Aj
        G_0
        tipo_modelo
    end
    
    methods
        %% ------------- Construtor para objetos fixos -------------- %%
        function obj = zener3D_fracionario_funcao_de_G_e_K_em_tensao(config)
            obj.k_bracos=config.k_bracos;
            obj.dimension=6;
            obj.flag_cinematica=1;
            obj.poisson=config.poisson;   
            obj.tipo_modelo='fracionario';
        end
        %% M�todo para fornecer vetor tempo 
        function value_n_time(obj,n_time)
                 obj.n_time=n_time;
                 obj.Aj=zeros(obj.n_time,obj.k_bracos);
        end 
        %% - M�todo para obter n�mero de variaveis internas e dimensao ---%%
        function size=get_size_var_inter(obj)
            size(1)=obj.dimension*obj.k_bracos;
            size(2)=obj.dimension;
        end
        
        %% ------ M�todo para atualizar variaveis de projeto ----------- %%
        function set_prop(obj,x)
            obj.G_inf=10^x(1);
            cont1=1;
            %% Par�metros viscoelasticos parcela desviadora
            obj.G_0=obj.G_inf;
            for i=1:obj.k_bracos
                obj.G(i)=10^x(cont1+i);
                obj.tau(i)=10^x(cont1+1+i);
                obj.alfa(i)=x(cont1+2+i);
                obj.G_0=obj.G_0+obj.G(i);
                cont1=cont1+2;
            end
             obj.K_inf=2*obj.G_0*(1+obj.poisson)/(3*(1-2*obj.poisson));
            
            %% Calculo Aj em posi��o invertida
            for i=1:obj.k_bracos
                obj.Aj(obj.n_time,i)=1; %inicializando o primeiro
                cont=obj.n_time-1;
                for k=1:obj.n_time-1
                    obj.Aj(cont,i)=obj.Aj(cont+1,i)*(k-1-obj.alfa(i))/k;
                    cont=cont-1;
                end
            end            
        end %Fim set prop
        
        
        %% ------------------ Fun��o modelo Zener 3D -------------------%%
        function[stress,var_inter,fail,et_1]=mat_solve(obj,result,it_atual)
           
            F_1=zeros(3);
            et_1=zeros(6,1);
            cont=1;
            for i=1:3
                for j=1:3
                    F_1(j,i)=result.F(cont,it_atual);
                    cont=cont+1;
                end
            end
            
            GradU_1=F_1-eye(3);
            def_infini_1=0.5*(GradU_1+GradU_1');
            for i=1:3
                et_1(i)=def_infini_1(i,i);
            end
            et_1(4)=def_infini_1(2,3); et_1(5)=def_infini_1(1,3);
            et_1(6)=def_infini_1(1,2);
            
            tn_1=result.t(it_atual);
            tn=result.t(it_atual-1);
            
            stress_iso_tot=zeros(6,1);
            sigma_inter=zeros(6,1);
            I=[1;1;1;0;0;0];
            dt=abs(tn_1-tn);
            %Contadores para operar com vetor hn e transformar em Hhn
            pa_ini=1;
            pa_fim=6;
            %Tra�o tensor deforma��o
            tr_e_1=(et_1(1)+et_1(2)+et_1(3));
            %Parcela volumetrica do tensor deforma��o
            vol_et_1=(1/3)*tr_e_1*I;
            %Parcela desviadora do tensor deforma��o
            dev_et_1=et_1-vol_et_1;
            %Parcela desviadora e volumetrica do tensor tens�o el�sico
            S_inf_vol=3*obj.K_inf*vol_et_1;
            S_inf_dev=2*obj.G_inf*dev_et_1;
            
            %% Calculo tens�o variaveis internas
            for j=1:obj.k_bracos
                %Calculo parametros
                B=1/(obj.tau(j)^obj.alfa(j));
                A=(dt^(-obj.alfa(j))*obj.tau(j)^obj.alfa(j)+1)*B;
                %Calculo derivada fracionaria
                Aj_1=obj.Aj((obj.n_time-it_atual+1):(obj.n_time-1),j);
                deri_frac_Q=GL_derivative3D(Aj_1,result.var_inter(pa_ini:pa_fim,1:(it_atual-1)),dt,obj.alfa(j));
                %Calculo valor variavel interna atual
                var_inter(pa_ini:pa_fim,1)=(1/A)*(B*2*obj.G(j)*dev_et_1-deri_frac_Q);
                %Atualizando vetores
                sigma_inter=sigma_inter+var_inter(pa_ini:pa_fim,1);
                pa_ini=pa_ini+6;  %posi�ao de inicio da variavel interna no vetor
                pa_fim=pa_fim+6;  %posi�ao final da variavel interna no vetor
                %% Valores elasticos desviadores de cada braco 
               stress_iso_tot=stress_iso_tot+2*obj.G(j)*dev_et_1;
                
            end                      
            stress=S_inf_vol+S_inf_dev+stress_iso_tot-sigma_inter;
            fail=false;
        end %Fim fun�ao do modelo constitutivo
        
        
    end %Fim fun�ao m�todos
end %Fim Classe