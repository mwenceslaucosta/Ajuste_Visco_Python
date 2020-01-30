function [delta_t_convergido,fail]=analise_conv_delta_t(config)
delta_t=[100 50 20 10 5 2];
n=length(delta_t);
fluencia.stress=4;
fluencia.tf=50E3;

mapeamento(1:500,1)=1:500;
v_convergencia=zeros(length(mapeamento),2);
eps=5E-5;
figure
for i=1:n
    v_tempo=0:delta_t(i):fluencia.tf;
    fluencia.deltat=delta_t(i);
    fluencia.t_carregamento=[1 2*60 fluencia.tf];
    config.fun_ensaio.set_ensaio_fluencia(fluencia);
    n_time=config.fun_ensaio.get_n_ensaio();
    config.fun_material.value_n_time(n_time)
    resultados_conv_1=Result();
    config.fun_material.set_prop(config.x0)
    resultados_conv_1.set_Result(config.fun_material,config.fun_ensaio);
    [resultados_conv_1,fail]=config.fun_ensaio.run(config.fun_material,resultados_conv_1);
    % Mapeamento para analise de convergencia
    for k=1:500
        for j=1:length(v_tempo)
            if v_tempo(j)==100*k
                mapeamento(k,2)=j;
                v_convergencia(k,2)=resultados_conv_1.e(1,j);
                break
            end
        end
    end
    erro=norm(v_convergencia(:,2)-v_convergencia(:,1));
    
    if erro<eps
        plot(v_tempo,resultados_conv_1.e(1,:))
        delta_t_convergido=delta_t(i-1);
        break
    end
    v_convergencia(:,1)=v_convergencia(:,2);
    plot(v_tempo,resultados_conv_1.e(1,:))
    hold on
end

legend('100 s','50 s','20 s','10 s','5 s')
%     switch n_models
%         case 1
%             title('Análise Convergência - Zener 3D Clássico ')
%         case 2
%              title('Análise Convergência - Zener 3D Fracionário')
%         case 3
%             title('Análise Convergência - Zener 3D Clássico Neo-Hookean')
%         case 4
%             title('Análise Convergência - Zener 3D Fracionário Neo-Hookean')
%     end

delta_t_convergido=min(delta_t_convergido);
end