classdef experimento_tracao < handle
   %% Classe para receber dados do ensaio de tração. 
   %% Recebe deformação na primeira coluna e tensão na primeira.
    properties
        t
        e
        stress
    end
    
    methods
        function set_experimento(obj,file_name)
            tabela=csvread(file_name);
            number_of_cols = size(tabela,2);
            if number_of_cols < 2
                error('colunas não confere')
            end
            obj.e=tabela(:,1);
            obj.stress=tabela(:,2);
            obj.t=tabela(:,1)/tabela(1,3); %Calcula tempo(def/taxa de def.)
          
        end
    end
    
end

