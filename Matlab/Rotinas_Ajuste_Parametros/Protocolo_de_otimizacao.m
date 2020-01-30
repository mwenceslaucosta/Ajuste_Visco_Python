function result=Protocolo_de_otimizacao(config,ensaios,experimentos)

%% Protocolo de otimizacao
%% Realiza 5 ajustes utilizando PSO
%% O melhor resultado da etapa anterior é utilizado em algoritmo baseado
%% em gradiente (BFGS).


%% Setup PSO

%1- Tamanho da população gerada a cada iteração
config_PSO.SwarmSize=150;
%2- Intervalo parâmetro de inercia w
config_PSO.InertiaRange=[0.5 0.5];
%3- Self adjustment parameter
config_PSO.SelfAdjustmentWeight=1.3;
%4- Social adjustment parameter
config_PSO.SocialAdjustmentWeight=1.3;
%5- Tolerancia criterio de parada
config_PSO.FunctionTolerance=2.5E-9;
%6- Numero maximo de iterações de Stall
config_PSO.MaxStallIterations=150;
%7 - Numero maximo de iterações
config_PSO.MaxIterations=150; 
%8- Numero de repetições do processo utilizando PSO
n_PSO=5;

%% SSE PSO
for i=1:n_PSO
    result.sse(i)=SSE(config.fun_material,ensaios,experimentos);
    
    %% Configuracão do problema de otimização
    funcao_objetivo{i}=@(x)result.sse(i).avaliar(x) ;
end

config_PSO.nVars=length(config.LB);
config_PSO.LB=config.LB;
config_PSO.UB=config.UB;

% Chamada PSO
for i=1:n_PSO
    display_1=['Iniciada ',num2str(i),'º iteração PSO'];
    disp(display_1)
    [x_PSO(i,:),f_PSO(i),hist_PSO{i}] = fun_runPSO(funcao_objetivo{i},config_PSO);
    
    %Salvando resultados em arquivo de texto
    iter=hist_PSO{i}.best(:,1);
    SSE_PSO=hist_PSO{i}.best(:,2);
    phi_PSO=hist_PSO{i}.phi_converg';
    T_PSO=table(iter,SSE_PSO,phi_PSO);
    writetable(T_PSO,(['SSE_PHI_PSO',num2str(i),'.csv']));
    
    %Salvando dados de tensao e deformacao do PSO para cada ensaio
    for j=1:length(ensaios)
        strain_ensaio(:,j)=result.sse(i).resultados(j).e(1,:);
        stress_ensaio(:,j)=result.sse(i).resultados(j).stress(1,:);
    end
    writetable(table(strain_ensaio),(['strain_PSO',num2str(i),'.csv']));
    writetable(table(stress_ensaio),(['stress_PSO',num2str(i),'.csv']));
end

%Salvando vetor tempo
data_time(:,1)=result.sse(1).resultados(1).t;
writetable(table(data_time),'time.csv');


[result.SSE_best_PSO,index]=min(f_PSO);
result.x_PSO=x_PSO;
result.f_PSO=f_PSO;
result.hist_PSO=hist_PSO;

result.best_x_PSO=x_PSO(index,:);
result.best_hist_PSO=hist_PSO{index};
result.x_otim=x_PSO(index,:);

%% Hibrido com o melhor resultado obtido pelo PSO.
result.BFGS.sse=SSE(config.fun_material,ensaios,experimentos);
funcao_objetivo_BFGS=@(x)result.BFGS.sse.avaliar(x) ;
x0=result.best_x_PSO;
LB=x0-1.5;
UB=x0+1.5;

if a=='fracionario'
    cont=0;
    for i=1:config.fun_material.k_bracos
    LB(cont+4)=0.01;
    UB(cont+4)=1;
    cont=cont+3;
    end 
end 
display_2='BFGS Iniciado';
disp(display_2)
options_BFGS = optimoptions('fmincon','Display','iter');
options_BFGS.OptimalityTolerance=1E-10;
[x_otim_BFGS,f_otim] =fmincon(funcao_objetivo_BFGS,x0,...
    [],[],[],[],LB,UB,[],options_BFGS);
result.x_otim=x_otim_BFGS;
result.SSE_otim=f_otim;

%Salvando dados de tensao e deformacao do BFGS
for j=1:length(ensaios)
    strain_ensaio(:,j)=result.BFGS.sse.resultados(j).e(1,:);
    stress_ensaio(:,j)=result.BFGS.sse.resultados(j).stress(1,:);
end

writetable(table(strain_ensaio),'strain_BFGS.csv');
writetable(table(stress_ensaio),'stress_BFGS.csv');

%Salvando vetor de propriedades
% Nomes_Colunas={'BFGS'};
Nomes_Colunas={'PSO1';'PSO2';'PSO3';'PSO4';'PSO5';'BFGS'};
T_xotim=table([x_PSO;x_otim_BFGS],'RowNames',Nomes_Colunas);
% T_xotim=table(x_otim_BFGS,'RowNames',Nomes_Colunas);

writetable(T_xotim,'x_otimos.csv');

end


