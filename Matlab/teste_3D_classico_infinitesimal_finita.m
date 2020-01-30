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
config.k_bracos=2;
config.poisson=0.46;
material_1=zener3D_ogden_Neo_Hookean(config);

delta_t=100;
[ensaios,experimentos] = Brant_fluencia_recu_UHMWPE(config,delta_t);

max_n_time=0; 
for i=1:length(ensaios)   
max_n_time=max(max_n_time,ensaios(i).get_n_ensaio());
end

x=[0.959474268362036,1.51569752353631,5.18653621382354,1.97476513530482,3.05140823297537];
result.sse1=SSE(material_1,ensaios,experimentos);
sse1=result.sse1.avaliar(x);

material_2=zener3D_classico_funcao_de_G_e_K_em_tensao(config);
result.sse2=SSE(material_2,ensaios,experimentos);
sse2=result.sse2.avaliar(x);

erro=abs(result.sse2.resultados(1, 1).e(1,:)-result.sse1.resultados(1, 1).e(1,:));
for i=1:3
p(1)=plot(result.sse1.resultados(i).t,result.sse1.resultados(i).e(1,:),'k');
hold on 
p(2)=plot(result.sse2.resultados(i).t,result.sse2.resultados(i).e(1,:),'k');
end 
