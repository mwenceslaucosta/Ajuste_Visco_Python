classdef experimento < handle
   
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
            obj.t=tabela(:,1);
            obj.e=tabela(:,2);
            if number_of_cols > 2 
                obj.stress=tabela(:,3); 
            end
        end
    end
    
end

