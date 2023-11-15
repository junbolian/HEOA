%_______________________________________________________________________________________%
% Human Evolutionary Optimization Algorithm (HEOA) source codes (version 1.0)           %
%                                                                                       %
%  Developed in MATLAB R2023a                                                           %
%  Author &  Inventor: Junbo Lian                                                       %
%  Corresponding author: Guohua Hui                                                     %
%  E-Mail:  deliver1982@163.com                                                         %
%_______________________________________________________________________________________%

function [avg_fitness_curve, Best_pos, Best_score, curve, search_history, fitness_history] = HEOA(N, Max_iter, lb, ub, dim, fobj)

jump_factor = abs(lb(1,1) - ub(1,1)) / 1000;

A = 0.6;  % Warning value
LN = 0.4; % Percentage of leaders
EN = 0.4; % Percentage of explorers
FN = 0.1; % Percentage of followers

% BestF: Best value in a certain iteration
% WorstF: Worst value in a certain iteration
% GBestF: Global best fitness value
% AveF: Average value in each iteration

LNNumber = round(N * LN); % Number of leaders
ENNumber = round(N * EN); % Number of explorers
FNNumber = round(N * FN); % Number of followers
if (max(size(ub)) == 1)
    ub = ub .* ones(1, dim);
    lb = lb .* ones(1, dim);
end

% Initialization
X0 = initializationLogistic(N, dim, ub, lb); %  Logistic chaotic mapping initialization representing the chaotic universe stage
X = X0;

% Compute initial fitness values
fitness = zeros(1, N);
for i = 1:N
    fitness(i) = fobj(X(i, :));
end
[fitness, index] = sort(fitness); % sort
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1); % Global best fitness value
AveF = mean(fitness);
for i = 1:N
    X(i, :) = X0(index(i), :);
end
curve = zeros(1, Max_iter);
avg_fitness_curve = zeros(1, Max_iter);
GBestX = X(1, :); % Global best position
X_new = X;
search_history = zeros(N, Max_iter, dim);
fitness_history = zeros(N, Max_iter);

% Start search
for i = 1:Max_iter
    disp(['Current number of iterations:', num2str(i)])
    BestF = fitness(1);
    WorstF = fitness(end);
    avg_fitness_curve(i) = AveF;
    R = rand(1);
    for j = 1:size(X, 1)
        % Human exploration stage-------------------------------------------------------------------------------------
        if i <= (1 / 4) * Max_iter
            X_new(j, :) = GBestX * (1 - i / Max_iter) + (mean(X(j, :)) - GBestX) * floor(rand() / jump_factor) * jump_factor + 0.2 * (1 - i / Max_iter) .* (X(j, :) - GBestX) .* Levy(dim); % Exploration stage
        else
            % Human development stage-------------------------------------------------------------------------------------
            for j = 1:LNNumber % Leaders
                if (R < A)
                    X_new(j, :) = 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * X(j, :) * exp((-i * randn(1)) / (rand(1) * Max_iter));
                else
                    X_new(j, :) = 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * X(j, :) + randn() * ones(1, dim);
                end
            end

            for j = LNNumber + 1:LNNumber + ENNumber % Explorers
                X_new(j, :) = randn() .* exp((X(end, :) - X(j, :)) / j^2);
            end

            for j = LNNumber + ENNumber + 1:LNNumber + ENNumber + FNNumber % Followers
                X_new(j, :) = X(j, :) + 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * rand(1, dim) .* (X(1, :) - X(j, :)); % Move with a randomly adaptive step size in the direction of the current best leader
            end

            for j = LNNumber + ENNumber + FNNumber + 1:N  % Losers
                X_new(j, :) = GBestX + (GBestX - X(j, :)) * randn(1);
            end
        end

        % Boundary control
        for j = 1:N
            for a = 1:dim
                if (X_new(j, a) > ub(a))
                    X_new(j, a) = ub(a);
                end
                if (X_new(j, a) < lb(a))
                    X_new(j, a) = lb(a);
                end
            end
        end

        % Update positions
        for j = 1:N
            fitness_new(j) = fobj(X_new(j, :));
        end
        for j = 1:N
            if (fitness_new(j) < GBestF)
                GBestF = fitness_new(j);
                GBestX = X_new(j, :);
            end
        end
        X = X_new;
        fitness = fitness_new;
        
        % Sorting and updating
        [fitness, index] = sort(fitness); % sort
        BestF = fitness(1);
        WorstF = fitness(end);
        for j = 1:N
            X(j, :) = X(index(j), :);
        end
        curve(i) = GBestF;
    end
    Best_pos = GBestX;
    Best_score = curve(end);
    search_history(:, i, :) = X;
    fitness_history(:, i) = fitness;
end


%  Levy spiral search strategy
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) *sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end


end
