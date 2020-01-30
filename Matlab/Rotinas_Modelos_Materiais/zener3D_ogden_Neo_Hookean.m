classdef zener3D_ogden_Neo_Hookean<modelos
    %% Modelo Viscoelástico 3D - cinemática finita 
    %% Ogden com N=1, recuperando modelo Neo-Hookean
    %% Utiliza modelo de Simo (para testar com parametros do Lucas Brant)
    
    %% ------------------- Informações modelo  --------------------- %%
    %k_bracos-> Numero de braços. Valor não atualizado durante otimização
    %x-> Vet. de projeto contendo parâmetros materiais.
    
    %x(1)=E_inf;  x(2)=E_1; x(3)=eta_1; x(4)=E2; x(5)=eta_2; ... ;
   
    
    properties
        Kvol
        E_inf
        E
        eta
        E0
        tau
        gamma
        k_bracos
        dimension
        n_time
        gamma_inf
        poisson
        G
        K     
        tipo_modelo
    end
    
    methods
        %% ------------- Construtor para objetos fixos -------------- %%
        function obj = zener3D_ogden_Neo_Hookean(config)
            obj.k_bracos=config.k_bracos;
            obj.dimension=6;
            obj.poisson=config.poisson;
            obj.tipo_modelo='classico';
        end
       %% Método para fornecer tamanho vetor tempo 
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
            obj.E_inf=10^x(1);
            obj.E0=obj.E_inf;
            cont1=1;
            for i=1:obj.k_bracos
                obj.E(i)=10^x(cont1+i);
                obj.eta(i)=10^x(cont1+i+1);
                obj.E0=obj.E0+obj.E(i);
                cont1=cont1+1;
            end
            obj.G=obj.E0/(2*(1+obj.poisson));
            obj.K=obj.E0/(3*(1-2*obj.poisson));
            %% Calculo gamma_i 
            for  i=1:obj.k_bracos
                obj.tau(i)=obj.eta(i)/obj.E(i);
                obj.gamma(i)=obj.E(i)/obj.E0;
            end
               obj.gamma_inf=obj.E_inf/obj.E0;
        end %Fim set prop
        
        
        %% --------------- Função modelo constitutivo -------------------%
        function[stress,var_inter,fail,et_1]=mat_solve(obj,result,it_atual)
           
            F_1=zeros(3);
            F_0=zeros(3);
            cont=1;
            for i=1:3
                for j=1:3
                    F_1(j,i)=result.F(cont,it_atual);
                    F_0(j,i)=result.F(cont,it_atual-1);
                    cont=cont+1;
                end
            end
            
            tn_1=result.t(it_atual);
            tn=result.t(it_atual-1);
            dt=abs(tn_1-tn); %Intervalo de tempo
            pa_ini=1;        %Contadores para operar com vetores de var.
            pa_fim=6;        %internas
            h_n_0=result.var_inter(:,it_atual-1);
 
            %% Calculo tensor C
            C=(F_1')*F_1;
            C_0=(F_0')*F_0;
            C_inv=inv(C);
            C_inv_0=inv(C_0);
            
            %% Tensor C e Cinv em vetor
            C_vector=zeros(6,1);
            C_vector_0=zeros(6,1);
            for i=1:3 
                C_vector(i,1)=C(i,i);
                C_vector_0(i,1)=C_0(i,i);
            end 
             C_vector(4,1)=C(2,3); C_vector_0(4,1)=C_0(2,3);
             C_vector(5,1)=C(1,3); C_vector_0(5,1)=C_0(1,3);
             C_vector(6,1)=C(1,2); C_vector_0(6,1)=C_0(1,2);
             
            Cinv_vec=zeros(6,1);
            Cinv_vec_0=zeros(6,1);
            for i=1:3 
                Cinv_vec(i,1)=C_inv(i,i);
                Cinv_vec_0(i,1)=C_inv_0(i,i);
            end 
            Cinv_vec(4,1)=C_inv(2,3); Cinv_vec_0(4,1)=C_inv_0(2,3);
            Cinv_vec(5,1)=C_inv(1,3); Cinv_vec_0(5,1)=C_inv_0(1,3);
            Cinv_vec(6,1)=C_inv(1,2); Cinv_vec_0(6,1)=C_inv_0(1,2);
     
            %Compatibilzando C para calculo correto da contração
            for j=4:6
                C_vector(j)=C_vector(j)*2;
                C_vector_0(j)=C_vector_0(j)*2;
            end
            
            %% Tensores desviadores
            I=[1;1;1;0;0;0];
            J=det(F_1);
            
            DEV_I=I-(1/3)*(dot(I,C_vector))*Cinv_vec;
            DEV_I_0=I-(1/3)*(dot(I,C_vector_0))*Cinv_vec_0;
            
            DEV_2_psi=obj.G*DEV_I;
            DEV_2_psi_0=obj.G*DEV_I_0;
            S_iso_n_1=DEV_2_psi;
            S_iso_n_0=DEV_2_psi_0;
            
            S_vol_n_1=obj.K*log(J)*Cinv_vec;
            S_hiperelastico=obj.gamma_inf*J^(-2/3)*S_iso_n_1+S_vol_n_1;
     
            sigma_var_inter=zeros(6,1);
            %% Calculo tensão variaveis internas
            for j=1:obj.k_bracos     
                %Calculo valor variavel interna atual
                H_n_0=h_n_0(pa_ini:pa_fim,1);
                H_n_1=exp(-dt/obj.tau(j))*H_n_0+exp(-dt/(2*obj.tau(j)))*(S_iso_n_1-S_iso_n_0);
                %Atualizando vetores
                var_inter(pa_ini:pa_fim,1)=H_n_1;
                pa_ini=pa_ini+6;  %posiçao de inicio da variavel interna no vetor
                pa_fim=pa_fim+6;  %posiçao final da variavel interna no vetor
                DEV_H=H_n_1-(1/3)*(dot(H_n_1,C_vector))*Cinv_vec;
                sigma_var_inter=sigma_var_inter+obj.gamma(j)*J^(-2/3)*DEV_H;      
            end
            
            
            %Tensão 2° Piola
            stress_S=S_hiperelastico+sigma_var_inter;
            %Retornando para tensorial
            stress_S_tensor(1,1)=stress_S(1); stress_S_tensor(2,2)=stress_S(2);
            stress_S_tensor(3,3)=stress_S(3); stress_S_tensor(2,3)=stress_S(4);
            stress_S_tensor(1,3)=stress_S(5); stress_S_tensor(1,2)=stress_S(6);
            stress_S_tensor(2,1)=stress_S_tensor(1,2);
            stress_S_tensor(3,2)=stress_S_tensor(2,3);
            stress_S_tensor(3,1)=stress_S_tensor(1,3);
            
            %Tensão Cauchy
            cauchy=(1/J)*F_1*stress_S_tensor*F_1';
            
            %Cauchy em vetor
            stress(1)=cauchy(1,1); stress(4)=cauchy(2,3);
            stress(2)=cauchy(2,2); stress(5)=cauchy(1,3);
            stress(3)=cauchy(3,3); stress(6)=cauchy(1,2);
            stress=stress';
            
            %Deformação verdadeira (logaritmica)
            b=F_1*F_1';
            V=sqrtm(b);
            def_ln=logm(V);
            
            et_1(1)=def_ln(1,1); et_1(4)=def_ln(2,3);
            et_1(2)=def_ln(2,2); et_1(5)=def_ln(1,3);
            et_1(3)=def_ln(3,3); et_1(6)=def_ln(1,2);       
            
            fail=false;
        end %Fim funçao do modelo constitutivo
    end %Fim métodos
end %Fim classe Neo-Hookean 3D