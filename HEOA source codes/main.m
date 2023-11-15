close all
clear
clc
N = 50; % Number of search agents
Function_name='F7'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration = 300; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj] = Get_Functions_details(Function_name);
[avg_fitness_curve, Best_pos, Best_score, Convergence1, search_history, fitness_history] = HEOA(N,Max_iteration,lb,ub,dim,fobj)

figure('Position',[454 445 1600 300])
%Draw search space
subplot(1,4,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])


%
subplot(1,4,2);
hold on
for k1 = 1:size(search_history, 1)
    for k2 = 1:size(search_history, 2)
        plot(search_history(k1, k2, 1), search_history(k1, k2, 2), '.', 'markersize', 1, 'MarkerEdgeColor', 'k', 'markerfacecolor', 'k');
    end
end
title('HEOA (x1 and x2 only)')
xlabel('x1')
ylabel('x2')
box on
axis tight




%
subplot(1,4,3);
hold on
plot(mean(fitness_history),'Linewidth', 1.5);
title('Average fitness')
xlabel('Iteration')
box on
axis tight



%Draw objective space
subplot(1,4,4);
semilogy(Convergence1,'Color','k','Linewidth', 1.5)
hold on
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');
axis tight
grid on
box on
legend('HEOA')



