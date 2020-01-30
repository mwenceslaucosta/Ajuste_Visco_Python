classdef Ensaios < handle
    
    properties
        mapeamento_experi %Tabela para relacionar experimento e numerico
    end
    
    methods (Abstract)
        [result,fail]=run(obj,mo_material,result)
        n=get_n_ensaio(obj) %Sendo fornecido no ensaio especifico
        sse=calc_sse(obj,experimento)
        time=get_time(obj)
        %fazer get t, para que nao seja necessario o ensaio ter exatamente
        %o objeto tempo. Desta forma deixo livre a forma que sera
        %forencido
        
    end
    
    methods
        
        function criar_mapeamento_experimentos(obj,time_experi)
            %Chamar
            %Verificacao do tempo de analise experimental e numerico
            time_numerico=obj.get_time();
            
            if max(time_experi)>max(time_numerico)%ficar atento max em matrizes
                error('Tempo de analise numerica menor que tempo experimental')
            end
            
            if min(time_experi)<min(time_numerico)%ficar atento max em matrizes
                error('Tempo de analise experimental menor que tempo numerio')
            end
            %% Cria vetor tempo corrigido, adicionando valores coincidentes com
            %% pontos experimentais
            %             ensaio.set_vetor_t_corrigido(time_experi)
%             for k=1:length(time_experi)
%                 for i=1:length(obj.t)-1
%                     if time_experi(k)>obj.t(i)&& time_experi(k)<obj.t(i+1)
%                         obj.t(i+1)= time_experi(k);
%                         obj.mapeamento_experi(k,1)=i+1;
%                         break
%                     end
%                 end
%             end
            % Metodo para criar tabela de mapemento entre valores experimentais
            % e valores numericos.
                        for i=1:length(time_experi)
                            for j=1:length( time_numerico)-1
                                if time_experi(i)>= time_numerico(j) && time_experi(i)<=time_numerico(j+1)
                                    obj.mapeamento_experi(i)=j;
                                    break %sai do loop menor, para nao continuar
                                end
                            end
                            %                 error('ponto nao encontrado')
                        end
            
        end %Fim metodo de mapeamento
    end %Fim metodos
end %Fim classe ensaios

