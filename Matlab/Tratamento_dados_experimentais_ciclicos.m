clear 
close
set_path 
dados_tensao=experimento();
dados_deformacao=experimento();
dados_tensao.set_experimento('stress_0_01%_s.csv')
dados_deformacao.set_experimento('strain_0_01%_s.csv')
time_tensao=dados_tensao.t;
time_deformacao=dados_deformacao.t;
mapeamento_deformacao=zeros(length(time_deformacao),1);

for i=1:length(time_tensao)
    for j=1:length( time_deformacao)-1
        if time_tensao(i)>= time_deformacao(j) && time_tensao(i)<=time_deformacao(j+1)
            mapeamento_deformacao(i)=j;
            break 
        end
    end
end

def_interpolado=zeros(length(time_deformacao),1);
for i=1:length(time_tensao)
    indice=mapeamento_deformacao(i);
    fator=(time_tensao(i)-time_deformacao(indice))/(time_deformacao(indice+1)-time_deformacao(indice));
    def_interpolado(i)=dados_deformacao.e(indice)+fator*(dados_deformacao.e(indice+1)-dados_deformacao.e(indice));
end

tabela=table(time_tensao,def_interpolado,dados_tensao.e);
writetable(tabela,'ciclico_0_1%_s.csv')
plot(def_interpolado,dados_tensao.e,'xb')
% t=27.817759599556567;
% t1=27.5940127029341;
% t2=27.8224816077991;
% e1=0.0160575837407595;
% e2=0.0144380136693702;
% m=(t-t1)/(t2-t1);
% e=e1+m*(e2-e1);