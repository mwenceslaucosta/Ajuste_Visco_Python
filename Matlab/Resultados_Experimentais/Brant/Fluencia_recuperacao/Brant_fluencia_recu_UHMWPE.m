function [ensaios,experimentos] = Brant_fluencia_recu_UHMWPE(config,delta_t)
%% Configura os ensaios para realizar ajustes das curvas experimentais 
% do Brant (somente parcela em fluencia)

%% Ensaios
vetor_ensaios=config.vetor_ensaios; 
n_ensaios=length(vetor_ensaios);
ensaios(n_ensaios)=config.fun_ensaio;
experimentos(n_ensaios)=experimento;
for i=1:n_ensaios
    switch (vetor_ensaios(i))
        case 1
%% Definindo ensaio 1
fluencia.stress=4;
fluencia.tf=165000.0;
fluencia.deltat=delta_t;        
fluencia.t_carregamento=[3 2*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('4MPa.csv')
%peso
% experimentos(i).e(1:36,2)=1;
% experimentos(i).e(27:30,2)=peso;
% experimentos(i).e(36:64,2)=peso;
         case 2
% Definindo ensaio 2
fluencia.stress=8;
fluencia.tf=165000.0;
fluencia.deltat=delta_t;
fluencia.t_carregamento=[3 4*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('Brant_8MPa.csv')
%Peso
% experimentos(i).e(1:43,2)=1;
% % experimentos(i).e(27:30,2)=peso;
% experimentos(i).e(43:70,2)=peso;

         case 3
% Definindo ensaio 3
fluencia.stress=16;
fluencia.tf=165000.0;
fluencia.deltat=delta_t;
fluencia.t_carregamento=[3 8*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('Brant_16MPa.csv')

    end %Fim Switch 
end %fim for 



end