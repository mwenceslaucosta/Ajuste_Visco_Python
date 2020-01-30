function deri_frac=GL_derivative3D(Aj_1,var_inter,dt,alfa)
%% Calculo derivada fracionária 3D
%Calculo por multiplicação de matrizes e sem redução de histórico
% Utiliza Aj calculados em posição invertida
deri_frac=dt^(-alfa)*var_inter*Aj_1;
end %Fim função para calculo da derivada fracionária 3D