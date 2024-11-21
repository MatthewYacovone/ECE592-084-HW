% A particle swarm optimizer
% to find the minimum/maximum of the MATLABs’ peaks function
% D---# of inputs to the function (dimension of problem)
clear; clc; close all;
%Parameters
ps=10;
D=2;
ps_lb=-3;
ps_ub=3;
vel_lb=-1;
vel_ub=1;
iteration_n = 50;
range = [-3, 3; -3, 3]; % Range of the input variables

upper = zeros(iteration_n, 1);
average = zeros(iteration_n, 1);
lower = zeros(iteration_n, 1);


% initialize population of particles and their velocities at time zero,
% format of pos= (particle#, dimension)
% construct random population positions bounded by VR

% need to bound positions
ps_pos=ps_lb + (ps_ub-ps_lb).*rand(ps,D);

% need to bound velocities between -mv,mv
ps_vel=vel_lb + (vel_ub-vel_lb).*rand(ps,D);

% initial pbest positions
p_best = ps_pos;

% returns column of cost values (1 for each particle)
f1='(ps_pos(i,2) - ps_pos(i,1)).^4 + 12*ps_pos(i,1)*ps_pos(i,2) - ps_pos(i,1) + ps_pos(i,2) - 3';
p_best_fit=zeros(ps,1);

for i=1:ps
    g1(i)=(ps_pos(i,2) - ps_pos(i,1)).^4 + 12*ps_pos(i,1)*ps_pos(i,2) - ps_pos(i,1) + ps_pos(i,2) - 3;
    p_best_fit(i)=g1(i);
end

p_best_fit;

hand_p3=plot3(ps_pos(:,1),ps_pos(:,2),p_best_fit','*k','markersize',15,'erase','xor');

% initial g_best
[g_best_val,g_best_idx] = max(p_best_fit);
%[g_best_val,g_best_idx] = min(p_best_fit); this is to minimize
g_best=ps_pos(g_best_idx,:);

% get new velocities, positions (this is the heart of the PSO algorithm)
for k=1:iteration_n
    for count=1:ps
    ps_vel(count,:) = 0.729*ps_vel(count,:)... % prev vel
    +1.494*rand*(p_best(count,:)-ps_pos(count,:))... % independent
    +1.494*rand*(g_best-ps_pos(count,:)); % social
    end
    ps_vel;

    % update new position
    ps_pos = ps_pos + ps_vel;

    %update p_best
    for i=1:ps
        g1(i)=(ps_pos(i,2) - ps_pos(i,1)).^4 + 12*ps_pos(i,1)*ps_pos(i,2) - ps_pos(i,1) + ps_pos(i,2) - 3;
        ps_current_fit(i)=g1(i);
        if ps_current_fit(i)>p_best_fit(i)
            p_best_fit(i)=ps_current_fit(i);
            p_best(i,:)=ps_pos(i,:);
        end
    end
    p_best_fit;

    %update g_best
    [g_best_val,g_best_idx] = max(p_best_fit);
    g_best=ps_pos(g_best_idx,:);
    % Fill objective function vectors
    upper(k) = max(p_best_fit);
    average(k) = mean(p_best_fit);
    lower(k) = min(p_best_fit);
end


figure(2);
x = 1:iteration_n;
plot(x, upper, 'o', x, average, 'x', x, lower, '*');
hold on;
plot(x, [upper average lower]);
hold off;
legend('Best', 'Average', 'Poorest');
xlabel('Iterations'); ylabel('Objective function value');