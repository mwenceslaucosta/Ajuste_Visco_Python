%% Rotina para ajuste de parâmetros experimentais do Lucas Brant em fluencia e recuperacao
clc
clear
close all

%% Adicionando caminhos de busca
set_path

%% Configuracao do problema para analise do intervalo de tempo

%Ensaio
config.fun_ensaio=ensaio_fluencia;
config.vetor_ensaios=[1 2 3];

%Modelo Material
config.k_bracos=1;
config.poisson=0.46;
config.fun_material=zener3D_classico_funcao_de_G_e_K_em_tensao(config);

%% Configuracão do problema de otimizacao
%Restricoes do problema de otimizacao
config.LB=[0 1 -4];
config.UB=[4.3 4.3 6];
config.x0=[log10(103.03) log10(143.70) 3.6045];

% Verificacao delta_t
% [delta_t_convergido,fail]=analise_conv_delta_t(config);
 delta_t_convergido=100;

%Definindo tempo máximo para ensaios numéricos 
[ensaios,experimentos] = Brant_fluencia_recu_UHMWPE(config,delta_t_convergido);
max_n_time=0; 
for i=1:length(ensaios)   
max_n_time=max(max_n_time,ensaios(i).get_n_ensaio());
end
n_time=max_n_time;
config.fun_material.value_n_time(n_time)
%Chamada protocolo de otimização
tic
result=Protocolo_de_otimizacao(config,ensaios,experimentos);
elapsedTime=toc;
%% pos-processamento
figure
for i=1:length(config.vetor_ensaios)
    if config.vetor_ensaios(i)==1
        p(1)=plot(result.BFGS.sse.resultados(i).t,result.BFGS.sse.resultados(i).e(1,:),'k');
    elseif config.vetor_ensaios(i)==2
        p(2)= plot(result.BFGS.sse.resultados(i).t,result.BFGS.sse.resultados(i).e(1,:),'--k');
    else
        p(3)= plot(result.BFGS.sse.resultados(i).t,result.BFGS.sse.resultados(i).e(1,:),'-.k');
    end
    hold on
    p(4)=plot(experimentos(i).t,experimentos(i).e,'bv');
    xlabel('tempo - [s]')
    ylabel('strain')
    xlim([0 165E3])
end

legend([p(1) p(2) p(3) p(4)],'4 MPa','8 MPa', '16 MPa','Experimental','Position',[330 220 0.005 0.005])
legend('boxoff')
disp(['Elapsed time - [s]:           ',num2str(elapsedTime)])
disp(['Melhor SSE PSO:               ',num2str(result.SSE_best_PSO)])
disp(['Melhores parametros PSO:      ',num2str(10.^result.best_x_PSO)])
disp(['SSE Hibrido:                  ',num2str(result.SSE_otim)])
disp(['Parametros Hibrido:           ',num2str(10.^result.x_otim)])



