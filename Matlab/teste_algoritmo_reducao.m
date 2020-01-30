clear
close
i_max=20;
k_max=20;
x=0:0.01:10;
delta_t=x(2)-x(1);
n_step=length(x);
alfa=0.5;
A_j=zeros(n_step,1);
A_j(1)=1;

cont=1;
for j=1:n_step-1
    %Calculo dos Aj até n_passos_toal
    A_j(j+1)=((j-1-alfa)/j)*A_j(cont);
    cont=cont+1;
end

%% Funcao para calculo dos fatores multiplicativos
%n_step: Numero de passos de tempo total do problema.

%vetor_1: Vetor contendando "k_max" valores da funcao teste.
vetor_1=ones(k_max,1);
%n_final: Numero de passos "n" apos o inicio da utilizacao da reducao
n_final=n_step-(i_max+k_max);
%Inicio da redução em (i_max+k_max+1)

%Chamada funcao para calculo dos fatores multiplicativos para obter Tn
fatores=zeros(k_max,n_final);
cont=1;
for n_passados=0:n_final
    fatores(:,cont)=calculo_fatores(n_passados,i_max,k_max,alfa);
    cont=cont+1;
end
%Calculo T0_1 para funcao teste f(t)=1.
T0_1=dot(vetor_1,fatores(:,1));

%Calculo T_inf_1 para funcao teste f(t)=1.
T_inf_1=sum(vetor_1);

%Calculo de Tn para funcao teste f(t)=1.
Tn_1=zeros(n_final+1,1);
for i=1:n_final+1
    Tn_1(i)=dot(vetor_1,fatores(:,i));
end

R_1=zeros(n_final+1,1);
%Calculo do coeficiente R_1(n) para funcao teste f(t)=1.
cont=1;
for i=1:n_final+1
    R_1(i)=(Tn_1(i)-T0_1)/(T_inf_1-T0_1);
end

%% Teste para funcao f(x)=x
%Utilizando todo historico
for i=1:n_step
    cont=i;
    summ=0;
    for j=1:i
        Aj_1=A_j(j);
        x_j=x(cont);
        summ=summ+Aj_1*x_j;
        cont=cont-1;
    end
    deri_total(i)=delta_t^(-alfa)*summ;
end

%% Derivada utilizando valores reduzidos
i_max_cont=i_max;
cont_n=1; %Utilizado para contabilizar steps apos inicio algoritmo de reducao
for i=1:n_step
    
    %% Ate i_max+k_max o calculo é feito pelo algoritmo clássico
    if i<=(i_max+k_max)
        cont=i;
        summ=0;
        for j=1:i
            Aj_1=A_j(j);
            x_j=x(cont);
            summ=summ+Aj_1*x_j;
            cont=cont-1;
        end
        deri_reduz(i)=delta_t^(-alfa)*summ;
    else
        
        %% Calculo historico mais recente
        
        %Calculo até i_max é obtido pelo algoritmo clássico
        summ=0;
        cont=i;
        for j=1:i_max
            Aj_1=A_j(j);
            x_j=x(cont);
            summ=summ+Aj_1*x_j;
            cont=cont-1;
        end
        deri_reduz_new=delta_t^(-alfa)*summ;
        
        %Calculo da influencia dos k valores anteriores apos o "i_max" considerado.
        fim_x_reduzido=cont;
        inico_x_reduzido=fim_x_reduzido-k_max+1;
        x_reduzido=x(inico_x_reduzido:fim_x_reduzido);
        [T0,T_inf]=get_T0_T_inf(x_reduzido,i_max,k_max,alfa);
        Aj_old=Aj_1(i_max+cont_n);
        Tn=(T0+R_1(cont_n)*(T_inf-T0));
        deri_reduz_old=delta_t^(-alfa)*Aj_old*Tn;  %Multiplicar pelo AJ!!!
        deri_reduz(i)=deri_reduz_old+deri_reduz_new; %
        cont_n=cont_n+1;
    end
    
end

x_reduzido=x(i_max:(i_max+k_max-1));

plot(x,deri_total)
%% Funcao para obter fatores multiplicativos utilizados no calculo de Tn
function fatores=calculo_fatores(n_passado,i_max,k_max,alfa)
cont=1;
n_posicoes=k_max;
fatores=zeros(n_posicoes,1);
fatores(1)=1;
for j=2:n_posicoes
    if j==2
        fatores(j)=(i_max+n_passado-alfa)/(i_max+n_passado+1);
    else
        fatores(j)=(i_max+n_passado-alfa+cont)/(i_max+n_passado+1+cont)*fatores(j-1);
        cont=cont+1;
    end
end
end

%% Funcao para obter T0 funcao a ser derivada

function [T0,T_inf]=get_T0_T_inf(vetor_func,i_max,k_max,alfa)
fatores_Tn_0_x(:,1)=calculo_fatores(0,i_max,k_max,alfa);

T0=dot(vetor_func,fatores_Tn_0_x);
T_inf=sum(vetor_func);
end

