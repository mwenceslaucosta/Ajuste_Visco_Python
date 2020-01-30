function deri_frac=GL_derivative3D(Aj_1,var_inter,dt,alfa)
%% Calculo derivada fracion�ria 3D
%Calculo por multiplica��o de matrizes e sem redu��o de hist�rico
% Utiliza Aj calculados em posi��o invertida
deri_frac=dt^(-alfa)*var_inter*Aj_1;
end %Fim fun��o para calculo da derivada fracion�ria 3D