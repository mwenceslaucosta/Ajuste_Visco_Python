classdef SSE < handle
       
    properties
        material
        ensaios
        experimentos
        resultados=Result.empty
        sse_ensaio
        sse0=1
    end
    
    methods
        
        %% Construtor
        function obj=SSE(material,ensaios,experimentos)
            
            if length(ensaios)~=length(experimentos)
                error('n. ensaios e experimentos não confere');
            end
            obj.material=material;
            obj.ensaios=ensaios;
            obj.experimentos=experimentos;
            
            n_ensaios=length(ensaios);
            obj.resultados(n_ensaios)=Result();
            obj.sse_ensaio=zeros(n_ensaios,1);
            
            for i=1:n_ensaios %Inicializando resultados para cada ensaio e criando mapeamento
                obj.resultados(i).set_Result(obj.material,obj.ensaios(i));
%               obj.ensaios(i).criar_mapeamento_experimentos(obj.experimentos(i).t)
                obj.ensaios(i).criar_mapeamento_experimentos(obj.experimentos(i).t)
            end
            
        end
        
        function set_normal(obj,sse0)
            obj.sse0=sse0;
        end
        
        %% Método para avaliar sse
        function sse=avaliar(obj,x)
            obj.material.set_prop(x);
            for i=1:length(obj.ensaios)
                [obj.resultados(i),fail]=obj.ensaios(i).run(obj.material,obj.resultados(i));
                if fail
                    error('ERRO DO ENSAIO!')
                end
                obj.sse_ensaio(i)=obj.ensaios(i).calc_sse(obj.experimentos(i),obj.resultados(i));
            end
            sse=sum(obj.sse_ensaio)/obj.sse0;
            
        end
    end
end

