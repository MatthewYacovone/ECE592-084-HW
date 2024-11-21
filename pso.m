% Particle swarm optimization
clear; close all; clc;

% Objective Function
function fitness_score = f(x1,x2)
fitness_score = (x2 - x1).^4 + 12*x1*x2 - x1 + x2 - 3;
end

% Parameters
popu = 10; % population size
D = 2; % num of inputs to the function (dimension of the problem)
iteration_n = 50; % num of iterations // termination criteria
w = 0.729;
c1 = 1.494;
c2 = 1.494;
pos_lb = -3;
pos_ub = 3;
vel_lb = -1;
vel_ub = 1;

% Initialize initial random positions and velocities
pos = pos_lb + (pos_ub-pos_lb).*rand(popu,D);
vel = vel_lb + (vel_ub-vel_lb).*rand(popu,D);

% Initialize containers
best = zeros(iteration_n, 1);
average = zeros(iteration_n, 1);
worst = zeros(iteration_n, 1);
p_best_fit=zeros(popu,1);

% Intialize personal bests
p_best = pos;

% Fitness scores of each initial personal best
for i=1:popu
    p_best_fit(i)=(pos(i,2) - pos(i,1)).^4 + 12*pos(i,1)*pos(i,2) - pos(i,1) + pos(i,2) - 3;
end

% Initialize global best
[g_best_fit,g_best_idx] = min(p_best_fit);
g_best=pos(g_best_idx,:);

% get new velocities, positions (this is the heart of the PSO algorithm)
for k=1:iteration_n
    for i=1:popu
        vel(i,:) = w*vel(i,:)... % prev vel
        +c1*rand*(p_best(i,:)-pos(i,:))... % independent
        +c2*rand*(g_best-pos(i,:)); % social
    end
    vel;

    % update new position
    pos = pos + vel;

    % update personal best
    for i=1:popu
        current_fit(i)=(pos(i,2) - pos(i,1)).^4 + 12*pos(i,1)*pos(i,2) - pos(i,1) + pos(i,2) - 3;
        if current_fit(i) < p_best_fit(i)
            p_best_fit(i)=current_fit(i);
            p_best(i,:)=pos(i,:);
        end
    end
    
    % updating global best
    [g_best_fit,g_best_idx] = min(p_best_fit);
    g_best=pos(g_best_idx,:);

    % updating best, average, and worst objective values
    best(k) = min(p_best_fit);
    average(k) = mean(p_best_fit);
    worst(k) = max(p_best_fit);
end

figure;
x = 1:iteration_n;
plot(x, best, 'ro', x, average, 'gx', x, worst, 'b*');
hold on;
plot(x, [best average worst]);
hold off;
legend('Best', 'Average', 'Worst');
xlabel('Generations'); ylabel('Objective Function Value'); title('Particle Swarm Optimization')







