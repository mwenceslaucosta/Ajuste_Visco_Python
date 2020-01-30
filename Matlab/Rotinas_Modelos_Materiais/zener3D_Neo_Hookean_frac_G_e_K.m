classdef zener3D_Neo_Hookean_frac_G_e_K<modelos
    %% Modelo Viscoelástico 3D - cinemática finita - Neo-Hookean 
    %% Formulado em termos de G e K.
    %% Parcela viscoelastica somente desviadora.
    %% ------------------- Informações modelo  --------------------- %%
    %k_bracos
    %x-> Vet. de projeto contendo parâmetros materiais.
    %x(1)=G_inf; x(2)=G_1; x(3)=tau1; x(4)=alfa_1; 
    %x(5)=G_2;   x(6)=tau2 ... x(7)=alfa_2 ;
 
    properties
        G
        G_inf
        K_inf
        tau
        alfa
        k_bracos
        dimension
        n_time
        Aj
        poisson        
        G_0
        tipo_modelo
    end
    
    methods
        %% ------------- Construtor para objetos fixos -------------- %%
        function obj = zener3D_Neo_Hookean_frac_G_e_K(config)
            obj.k_bracos=config.k_bracos;
            obj.dimension=6;
            obj.poisson=config.poisson;
             obj.tipo_modelo='fracionario';
        end
        
        %% Método para fornecer tamanho vetor tempo 
        function value_n_time(obj,n_time)
                 obj.n_time=n_time;
                 obj.Aj=zeros(obj.n_time,obj.k_bracos);
        end
        
        %% - Método para obter número de variaveis internas e dimensao ---%%
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
                obj.alfa(i)=x(cont1+i+2);
                obj.G_0=obj.G_0+obj.G(i);
                cont1=cont1+2;
            end
          obj.K_inf=2*obj.G_0*(1+obj.poisson)/(3*(1-2*obj.poisson));
         
        
        %% Calculo Aj em posição invertida
            for  i=1:obj.k_bracos
                obj.Aj(obj.n_time,i)=1; %inicializando primeiro
                cont=obj.n_time-1;
                for k=1:obj.n_time-1
                    obj.Aj(cont,i)=obj.Aj(cont+1,i)*(k-1-obj.alfa(i))/k;
                    cont=cont-1;
                end
            end
        end %Fim set prop
        
        
        %% --------------- Função modelo constitutivo -------------------%
        function[stress,var_inter,fail,et_1]=mat_solve(obj,result,it_atual)
 
            %% Realocando tensor F para obter C em notação tensorial
            %F é recebio em notação vetorial de 9 componentes conforme JOG
            F_1=zeros(3);
            cont=1;
            for i=1:3
                for j=1:3
                    F_1(j,i)=result.F(cont,it_atual);
                    cont=cont+1;
                end
            end
            
            tn_1=result.t(it_atual);
            tn=result.t(it_atual-1);
            dt=abs(tn_1-tn); %Intervalo de tempo
            pa_ini=1;        %Contadores para operar com vetores de var.
            pa_fim=6;        %internas
            Q_tot=zeros(6,1);
            S_iso_tot=zeros(6,1);
            
            %% Calculo tensor C
            C=(F_1')*F_1;
            C_inv=inv(C);
            
            %% Tensor C e Cinv em vetor
            C_vector=zeros(6,1);
            for i=1:3 
                C_vector(i,1)=C(i,i);
            end 
             C_vector(4,1)=C(2,3);
             C_vector(5,1)=C(1,3);
             C_vector(6,1)=C(1,2);
             
            Cinv_vec=zeros(6,1);
            for i=1:3 
                Cinv_vec(i,1)=C_inv(i,i);
            end 
            Cinv_vec(4,1)=C_inv(2,3);
            Cinv_vec(5,1)=C_inv(1,3);
            Cinv_vec(6,1)=C_inv(1,2);
     
            %Compatibilzando C para calculo correto da contração
            for j=4:6
                C_vector(j)=C_vector(j)*2;
            end
            
            %% Tensores desviadores
            I=[1;1;1;0;0;0];
            J=det(F_1);
        
            DEV_I=I-(1/3)*(dot(I,C_vector))*Cinv_vec;
            %DEV_Ciso=Ciso-(1/3)*(dot(Ciso,C_voigt))*Cinv_voigt;
            
            DEV_2_psi_inf=DEV_I*obj.G_inf;
            S_iso_inf=J^(-2/3)*DEV_2_psi_inf;
            S_vol_inf=obj.K_inf*log(J)*Cinv_vec;
            S_hiperelastico_inf=S_iso_inf+S_vol_inf;            
           
            for j=1:obj.k_bracos
                %% Calculo tensão variaveis internas
                %Calculo parametros
                 B=1/(obj.tau(j)^obj.alfa(j));
                 A=(dt^(-obj.alfa(j))*obj.tau(j)^obj.alfa(j)+1)*B;
                
                %Calculo derivada fracionaria
                Aj_1=obj.Aj((obj.n_time-it_atual+1):(obj.n_time-1),j);
                deri_frac_Q=GL_derivative3D(Aj_1,result.var_inter(pa_ini:pa_fim,1:(it_atual-1)),dt,obj.alfa(j));
                %Calculo valor variavel interna atual
                DEV_2_psi_i=DEV_I*obj.G(j)*J^(-2/3);
                Q_n_1=(1/A)*(B*DEV_2_psi_i-deri_frac_Q);
                
                %Atualizando vetores
                var_inter(pa_ini:pa_fim,1)=Q_n_1;
                pa_ini=pa_ini+6;  %posiçao de inicio da variavel interna no vetor
                pa_fim=pa_fim+6;  %posiçao final da variavel interna no vetor
                Q_tot=Q_n_1+Q_tot;
                
                %% Valores elasticos desviadores de cada braco 
                S_iso_tot=S_iso_tot+DEV_2_psi_i;
                
            end
                        
            %Tensão 2° Piola
            stress_S=S_hiperelastico_inf+S_iso_tot-Q_tot;
            %Retornando para tensorial
            %stress_S=[S_11 S_22 S_33 S_23 S_13 S_11]
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
            
            %Piola
            %Tensão Cauchy
            cauchy=(1/J)*F_1*stress_S_tensor*F_1';
            
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
            
            fail=false;
        end %Fim funçao do modelo constitutivo
    end %Fim métodos
end %Fim classe Neo-Hookean 3D