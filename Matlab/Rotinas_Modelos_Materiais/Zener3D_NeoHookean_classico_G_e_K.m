classdef Zener3D_NeoHookean_classico_G_e_K<modelos
    %% Rotina modelo Material de Zener clássico, em cinematica finita, uti-
    %% lizando funcao de energia de Neo Hookean. Formulado em termos de G e K
    %% para cada braço.
    %% Parcela viscoelastica somente desviadora. 
    % x-> Vet. de projeto contendo parâmetros materiais.
    % x(1)=G_inf; x(2)=G_1; x(3)=tau1; x(4)=G_2; x(5)=tau2 ... 
    
    properties
        G
        K
        G_inf
        K_inf
        tau
        k_bracos
        dimension
        n_time
        poisson
        G_0
        tipo_modelo
    end
    
    methods
        %% ------------- Construtor para objetos fixos -------------- %%
        function obj = Zener3D_NeoHookean_classico_G_e_K(config)
            obj.k_bracos=config.k_bracos;
            obj.dimension=6;
            obj.poisson=config.poisson;
            obj.tipo_modelo='classico';
        end
        
        %% Método para fornecer tamanho vetor tempo 
        function value_n_time(obj,n_time)
                 obj.n_time=n_time;
        end
        
        %% - Método para obter número de variaveis internas e dimensao -%%
        function size=get_size_var_inter(obj)
            size(1)=obj.dimension*obj.k_bracos;
            size(2)=obj.dimension;
        end
        %% ------ Método para atualizar variaveis de projeto ----------- %%
        function set_prop(obj,x)
            obj.G_inf=10^x(1);
            obj.G_0=obj.G_inf;
            cont1=1;
            for i=1:obj.k_bracos
                obj.G(i)=10^x(cont1+i);
                obj.tau(i)=10^x(cont1+i+1);
                obj.G_0=obj.G_0+obj.G(i);
                cont1=cont1+1;
            end
            obj.K_inf=2*obj.G_0*(1+obj.poisson)/(3*(1-2*obj.poisson));
        end %Fim set prop
        
   %% --------------- Função modelo constitutivo -------------------%
        function[stress,var_inter,fail,et_1]=mat_solve(obj,result,it_atual)
 
            %F é recebio em notação vetorial de 9 componentes conforme JOG
            F_1=zeros(3);
            F_0=zeros(3);
            cont=1;
            tn_1=result.t(it_atual);
            tn=result.t(it_atual-1);
            hn_0=result.var_inter(:,it_atual-1);
            dt=abs(tn_1-tn); %Intervalo de tempo
            
            %% Realocando tensor F para obter C em notação tensorial
            %F=[F11 F21 F31 F12 F22 F32 F13 F23 F33];
            for i=1:3
                for j=1:3
                    F_1(j,i)=result.F(cont,it_atual);
                    F_0(j,i)=result.F(cont,it_atual-1);
                    cont=cont+1;
                end
            end
            Grad_U=F_1-eye(3);
            def_inf=0.5*(Grad_U'+Grad_U);
            pa_ini=1;        %Contadores para operar com vetores de var.
            pa_fim=6;        %internas
            
            %% Calculo tensor C
            C_1=(F_1')*F_1;
            C_inv_1=inv(C_1);
           
            C_0=(F_0')*F_0;
            C_inv_0=inv(C_0);
            
            %% Tensor C e Cinv em vetor (C11 C22 C33 C23 C13 C12)
            C_vector_1=zeros(6,1);
            C_vector_0=zeros(6,1);
            
            for i=1:3 
                C_vector_1(i,1)=C_1(i,i);
                C_vector_0(i,1)=C_0(i,i);
                def_inf_vec(i,1)=def_inf(i,i);
            end
            
            C_vector_1(4,1)=C_1(2,3); C_vector_0(4,1)=C_0(2,3);
            C_vector_1(5,1)=C_1(1,3); C_vector_0(5,1)=C_0(1,3);
            C_vector_1(6,1)=C_1(1,2); C_vector_0(6,1)=C_0(1,2);
            
            def_inf_vec(4,1)=def_inf(2,3);
            def_inf_vec(5,1)=def_inf(1,3);
            def_inf_vec(6,1)=def_inf(1,2);
            
            Cinv_vec_1=zeros(6,1);
            Cinv_vec_0=zeros(6,1);
            
            for i=1:3 
                Cinv_vec_1(i,1)=C_inv_1(i,i);
                Cinv_vec_0(i,1)=C_inv_0(i,i);
            end 
            Cinv_vec_1(4,1)=C_inv_1(2,3); Cinv_vec_0(4,1)=C_inv_0(2,3);
            Cinv_vec_1(5,1)=C_inv_1(1,3); Cinv_vec_0(5,1)=C_inv_0(1,3);
            Cinv_vec_1(6,1)=C_inv_1(1,2); Cinv_vec_0(6,1)=C_inv_0(1,2);
     
            %Compatibilzando C para calculo correto da contração
            for j=4:6
                C_vector_1(j)=C_vector_1(j)*2;
                C_vector_0(j)=C_vector_0(j)*2;
            end
            
            %% Tensores desviadores
            sigma_inter_dev=zeros(6,1);
            I=[1;1;1;0;0;0];
            J_1=det(F_1);
            J_0=det(F_0);
        
            DEV_I_1=I-(1/3)*(dot(I,C_vector_1))*Cinv_vec_1;
            DEV_I_0=I-(1/3)*(dot(I,C_vector_0))*Cinv_vec_0;
            
            S_iso_inf=J_1^(-2/3)*DEV_I_1*obj.G_inf;
            S_vol_inf=obj.K_inf*J_1*(J_1-1)*Cinv_vec_1;
%           S_vol_inf=obj.K_inf*log(J_1)*Cinv_vec_1;
            S_hiperelastico_inf=S_iso_inf+S_vol_inf;            
           
            for j=1:obj.k_bracos
   %% Calculo tensão variaveis internas -  Parcela viscoelástica desviadora
                Hhn_0=hn_0(pa_ini:pa_fim,1);
                A=exp(-dt/obj.tau(j))*Hhn_0;
                S_1=J_1^(-2/3)*obj.G(j)*DEV_I_1;
                S_0=J_0^(-2/3)*obj.G(j)*DEV_I_0;
                B=exp(-dt/(2*obj.tau(j)))*(S_1-S_0);
                Hhn_1=A+B;
                sigma_inter_dev=sigma_inter_dev+Hhn_1;
                var_inter(pa_ini:pa_fim,1)=Hhn_1;
                pa_ini=pa_ini+6;  %posiçao inicial no vetor var_inter
                pa_fim=pa_fim+6;  %posiçao final no vetor var_inter
                
            end             
            %Tensão 2° Piola
            stress_S=S_hiperelastico_inf+sigma_inter_dev;
            
            %Retornando para tensorial
            %stress_S=[S_11 S_22 S_33 S_23 S_13 S_12]
            stress_S_tensor=zeros(3);
            for i=1:3 
                stress_S_tensor(i,i)=stress_S(i);
            end 
            stress_S_tensor(1,2)=stress_S(6);
            stress_S_tensor(2,3)=stress_S(4);
            stress_S_tensor(1,3)=stress_S(5); 
            
            stress_S_tensor(2,1)=stress_S_tensor(1,2);
            stress_S_tensor(3,2)=stress_S_tensor(2,3);
            stress_S_tensor(3,1)=stress_S_tensor(1,3);
            
            %Tensão Cauchy
            cauchy=(1/J_1)*F_1*stress_S_tensor*F_1';
            
            %Cauchy em vetor
            stress=zeros(6,1);
            for i=1:3 
                stress(i,1)=cauchy(i,i);
            end 
            stress(4)=cauchy(2,3);
            stress(5)=cauchy(1,3);
            stress(6)=cauchy(1,2);
        
            %Deformação verdadeira (logaritmica)
            b=F_1*F_1';
            V=sqrtm(b);
            def_ln=logm(V);
            
            et_1=zeros(6,1);
            for i=1:3 
                et_1(i,1)=def_ln(i,i);
            end 
            
            et_1(4)=def_ln(2,3);
            et_1(5)=def_ln(1,3);
            et_1(6)=def_ln(1,2);
            
            et_1=def_inf_vec;
            fail=false;
        end %Fim funçao do modelo constitutivo
    end %Fim métodos
end %Fim classe 

