classdef modelos<handle
    
    properties
        
    end
    
    methods
        
    end
    methods (Abstract)
        %valores que todos os modelos devem fazer
        %caso contrario - error
        [stress,var_inter,fail]=mat_solve(obj,result,i)
        %result; time, stress, strain
        
        %O ensaio chama modelos e fornece tais entradas
        size_var_inter=get_size_var_inter(obj)
        set_prop(obj,x)
        
    end
    
end

