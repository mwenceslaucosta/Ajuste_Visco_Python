function [x_PSO,f_PSO,hist_PSO] = fun_runPSO(funcao_objetivo,config_PSO)
hist_PSO = [];
cont=1;
delta_best_worst_0=0;
eps=config_PSO.FunctionTolerance;
ns=0.9; %Amostragem para verificar max e min (90%)

fig1=figure;
fig2=figure;
%% Chamada Otimizador

options_swarm = optimoptions('particleswarm','Display','iter','OutputFcn',@nested_PSO_hist);
options_swarm.SwarmSize=config_PSO.SwarmSize;
options_swarm.InertiaRange=config_PSO.InertiaRange;
options_swarm.SelfAdjustmentWeight=config_PSO.SelfAdjustmentWeight;
options_swarm.SocialAdjustmentWeight=config_PSO.SocialAdjustmentWeight;
options_swarm.MaxIterations=config_PSO.MaxIterations;
options_swarm.UseParallel=true;	
options_swarm.MaxStallIterations=config_PSO.MaxStallIterations;
LB=config_PSO.LB;
UB=config_PSO.UB;
n_values=ceil(options_swarm.SwarmSize*ns);

[x_PSO, f_PSO]=particleswarm(funcao_objetivo,config_PSO.nVars,LB,UB,options_swarm);

%% Nested function para conseguir salvar historico do algoritmo e plotar
%% outras formas de graficos
    function stop = nested_PSO_hist(optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                hist_PSO.best(cont,1)=cont-1;
                hist_PSO.best(cont,2)=optimValues.bestfval;
                
                %% Plot iter x SSE
                set(groot,'CurrentFigure',fig1);
                semilogy(cont-1,optimValues.bestfval,'ko')
                hold on
                xlabel('Iteração')
                ylabel('SSE')
                drawnow;

                swarm_crescente_ordem=sort(optimValues.swarmfvals);
                swarm_crescente_n_values=swarm_crescente_ordem(1:n_values);
                best_0=min(swarm_crescente_n_values);
                worst_0=max(swarm_crescente_n_values);
                delta_best_worst_0=worst_0-best_0;
                
                %% Plot iter x phi_convergencia
                set(groot,'CurrentFigure',fig2);
                delta_best_worst_1=phi_conv_1(optimValues);
                phi_converg=delta_best_worst_1/delta_best_worst_0;
                semilogy(cont-1,phi_converg,'ko')
                hold on
                xlabel('Iteração')
                ylabel('Phi_k')
                drawnow;
                
                hist_PSO.phi_converg(cont)=phi_converg;
                cont=cont+1;
            case 'iter'
                hist_PSO.best(cont,1)=cont-1;
                hist_PSO.best(cont,2)=optimValues.bestfval;
                
                vec_optimValues=[hist_PSO.best(cont-1,2) optimValues.bestfval];
                set(groot,'CurrentFigure',fig1);
                semilogy(cont-1,optimValues.bestfval,'ko')
                semilogy([cont-2 cont-1],vec_optimValues,'k')
                xlim([0 cont-1])
                drawnow;
                
                %% Verificacao convergecia Swarm (diferença relativa)
                delta_best_worst_1=phi_conv_1(optimValues);
                phi_converg=delta_best_worst_1/delta_best_worst_0;
                vec_phi_converg=[hist_PSO.phi_converg(cont-1) phi_converg];
                set(groot,'CurrentFigure',fig2);
                semilogy(cont-1,phi_converg,'ko')
                semilogy([cont-2 cont-1],vec_phi_converg,'k')
                xlim([0 cont-1])
                drawnow;
                if (cont-1)>100
                if phi_converg<=eps
                    stop=true;
                    disp(['PSO convergiu. Phi_conv=', num2str(phi_converg)])
                end
                end
                           
                hist_PSO.phi_converg(cont)=phi_converg;
                cont=cont+1;
                
            case 'done'
                hold off
            otherwise
        end
    end

    function delta_best_worst_1=phi_conv_1(optimValues)
        swarm_crescente_ordem_1=sort(optimValues.swarmfvals);
        swarm_crescente_n_values_1=swarm_crescente_ordem_1(1:n_values);
        best_1=min(swarm_crescente_n_values_1);
        worst_1=max(swarm_crescente_n_values_1);
        delta_best_worst_1=worst_1-best_1;
    end

end