function resultado_protocolo = teste_protocolo(config)

%% Etapa 1 - Criação de resultados artificiais
% Vamos criar um conjunto de resultados experimentais fictícios com o mesmo
% número de ensaios fornecidos
ensaios=config.ensaios;
n_ensaios=length(ensaios);
experimentos(n_ensaios)=experimento();

% Atribuir as propriedades de referencia ao material

material=config.fun_material;
propriedades_referencia=config.x0;
material.set_prop(propriedades_referencia)

for iexp = 1:n_ensaios
    % Criando uma estrutura para guardar os resultados temporários
    % e configurando com o material e ensaio de interesse
    resultados = Result();
    resultados.set_Result(material, ensaios(iexp))
    
    % Rodar o ensaio
    [resultados, falhou] = ensaios(iexp).run(material, resultados);
    if falhou
        error(['Criação de resultado falhou. id=', num2str(iexp)])
    end
    
    % Com os resultados, vamos "enganar" o protocolo criando um experimento
    experimentos(iexp).t = resultados.t;
    experimentos(iexp).e = resultados.e(1,:)';
    experimentos(iexp).stress = resultados.stress(1,:)';
end

%% Etapa 2 - Utilização do Protocolo de ajuste para recuperar o material
% A filosofia desta etapa é testar o protocolo de ajuste.
% Vamos solicitar ao protocolo que ele tente ajustar os nossos resultados
% experimentais fictícios. Com este procedimento esperamos que o protocolo
% encontre as mesmas propriedades de referencia (propriedades_referencia)
% que utilizamos para gerar os nossos "experimentos"

% Chamada do Protocolo
resultado_protocolo = Protocolo_de_otimizacao(config, ensaios, experimentos);

%% Etapa 3 - Verificação da resposta do protocolo
% Vamos utilizar todos os resultados fornecidos pelo protocolo e comparar
% com a solução que nós conhecemos e esperamos. Sabemos que o conjunto de
% parâmetros ótimo é o próprio vetor "propriedades_referencia", e que o SSE
% esperado é zero, isto é, temos capacidade de obter a solução EXATA do
% problema.
disp([' SSE: ',num2str(resultado_protocolo.SSE_otim)]) 
propriedades_referencia=config.x0;
disp([' Parâmetros encontrados:    ', ...
    num2str(resultado_protocolo.best_x_PSO)])
disp([' Parâmetros de referência:  ', ...
    num2str(propriedades_referencia)])
dif=abs(resultado_protocolo.x_otim-propriedades_referencia);
disp([' |Ref-obtidos|:             ', num2str(dif)])

end