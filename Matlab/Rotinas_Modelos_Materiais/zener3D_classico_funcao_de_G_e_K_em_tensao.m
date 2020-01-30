classdef zener3D_classico_funcao_de_G_e_K_em_tensao<modelos
    
    %% --------------------Informações modelo zener 3D --------------------- %%
    %Modelo formulado com funções energia de deformação separada para
    %cada braço, em termos de G e K.
    %parcela viscoelastica somente desviadora. 
    %Modelo criado para verificação da formulação em tensão, para ficar
    %consistente com os outros modelos 3D apresentados na dissertação
    properties
        G_inf
        K_inf
        K_0
        G_0
        G
        K
        tau_dev        
        k_bracos
        dimension
        flag_cinematica
        n_time
        Algoritmo
        poisson
        tipo_modelo           
    end
    
    methods
        %% ------------- Construtor para objetos fixos -------------- %%
        function obj = zener3D_classico_funcao_de_G_e_K_em_tensao(config)
            obj.k_bracos=config.k_bracos;
            obj.dimension=6;
            obj.flag_cinematica=1;   
            obj.poisson=config.poisson;
            obj.tipo_modelo='classico';
        end
             %% Método para fornecer vetor tempo 
        function value_n_time(obj,n_time)
                 obj.n_time=n_time;
        end    
        %% - Método para obter número de variaveis internas e dimensao ---%%
        function size=get_size_var_inter(obj)
            size(1)=obj.dimension*obj.k_bracos;
            size(2)=obj.dimension;
        end
        
        %% ------ Método para atualizar variaveis de projeto ----------- %%
        function set_prop(obj,x)
            
            obj.G_inf=10^x(1);
            cont1=1;            
            %% Parâmetros viscoelasticos parcela desviadora
            obj.G_0=obj.G_inf;
            for i=1:obj.k_bracos
                obj.G(i)=10^x(cont1+i);
                obj.tau_dev(i)=10^x(cont1+1+i);
                obj.G_0=obj.G_0+obj.G(i);
                cont1=cont1+1;
            end
            %parametro volumetrico para hipotese de parcela
            %elastica puramente volumetrica
            obj.K_inf=2*obj.G_0*(1+obj.poisson)/(3*(1-2*obj.poisson));
            
        end %Fim set prop
        
        
        %% ------------------ Função modelo Zener 3D -------------------%%
        function[stress,var_inter,fail,et_1]=mat_solve(obj,result,it_atual)
            F_1=zeros(3);
            F_0=zeros(3);
            et_1=zeros(6,1);
            et_0=zeros(6,1);
            tn_1=result.t(it_atual);
            tn=result.t(it_atual-1);
            hn_0=result.var_inter(:,it_atual-1);
            cont=1;
          
            for i=1:3
                for j=1:3
                    F_1(j,i)=result.F(cont,it_atual);
                    F_0(j,i)=result.F(cont,it_atual-1);
                    cont=cont+1;
                end
            end
            
            GradU_1=F_1-eye(3);
            GradU_0=F_0-eye(3);
            def_infini_1=0.5*(GradU_1+GradU_1');
            def_infini_0=0.5*(GradU_0+GradU_0');
            
            for i=1:3
                et_1(i)=def_infini_1(i,i);
                et_0(i)=def_infini_0(i,i);
            end
            et_1(4)=def_infini_1(2,3); et_1(5)=def_infini_1(1,3);
            et_1(6)=def_infini_1(1,2); et_0(4)=def_infini_0(2,3);
            et_0(5)=def_infini_0(1,3); et_0(6)=def_infini_0(1,2);
            
            sigma_inter_dev=zeros(6,1);
            I=[1;1;1;0;0;0];
            dt=abs(tn_1-tn);
            
            %Contadores para operar com vetor hn e transformar em Hhn
            pa_ini=1;
            pa_fim=6;
            %Traço tensor deformação
            tr_e_1=(et_1(1)+et_1(2)+et_1(3));
            tr_e_0=(et_0(1)+et_0(2)+et_0(3));
            %Parcela volumetrica do tensor deformação
            vol_et_1=(1/3)*tr_e_1*I;
            vol_et_0=(1/3)*tr_e_0*I;
            %Parcela desviadora do tensor deformação
            dev_et_1=et_1-vol_et_1;
            dev_et_0=et_0-vol_et_0;
            %Parcela desviadora e volumetrica do tensor tensão elásico
            S_inf_vol=3*obj.K_inf*vol_et_1;
            S_inf_dev=2*obj.G_inf*dev_et_1;
            
            %Parcela viscoelástica desviadora
            for j=1:obj.k_bracos
                Hhn_0=hn_0(pa_ini:pa_fim,1);
                A=exp(-dt/obj.tau_dev(j))*Hhn_0;
                B=exp(-dt/(2*obj.tau_dev(j)))*2*obj.G(j)*(dev_et_1-dev_et_0);
                Hhn_1=A+B;
                sigma_inter_dev=sigma_inter_dev+Hhn_1;
                var_inter(pa_ini:pa_fim,1)=Hhn_1;
                pa_ini=pa_ini+6;  %posiçao inicial no vetor var_inter
                pa_fim=pa_fim+6;  %posiçao final no vetor var_inter
            end
            
            %Tensão total
            stress=S_inf_vol+S_inf_dev+sigma_inter_dev;
            fail=false;
        end %Fim funçao do modelo constitutivo
    end
end