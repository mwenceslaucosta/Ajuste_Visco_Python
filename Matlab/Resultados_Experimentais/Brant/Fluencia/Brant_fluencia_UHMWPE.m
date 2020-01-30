function [ensaios,experimentos] = Brant_fluencia_UHMWPE(config,delta_t)
%% Configura os ensaios para realizar ajustes das curvas experimentais 
% do Brant

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
fluencia.tf=86E3;
fluencia.deltat=delta_t;        
fluencia.t_carregamento=[1 2*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('Brant_4MPa_creep.csv')
         case 2
% Definindo ensaio 2
fluencia.stress=8;
fluencia.tf=86E3;
fluencia.deltat=delta_t;
fluencia.t_carregamento=[1 4*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('Brant_8MPa_creep.csv')
         case 3
% Definindo ensaio 3
fluencia.stress=16;
fluencia.tf=86E3;
fluencia.deltat=delta_t;
fluencia.t_carregamento=[1 8*60 24*60*60];
ensaios(i).set_ensaio_fluencia(fluencia);
experimentos(i).set_experimento('Brant_16MPa_creep.csv')
    end %Fim Switch 
end %fim for 



end