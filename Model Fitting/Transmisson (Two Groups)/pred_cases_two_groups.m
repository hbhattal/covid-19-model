function matrix = pred_cases_two_groups(T, D, N, para, tspan, y1)

% >>> define beta matrix in terms of para

% no mixing, varying contact rates between groups
% beta = [para(1), 0; 0, para(2)];

% homogenous mixing, varying contact rates between groups
% beta = [para(1)*(N(1)/sum(N)), para(1)*(N(2)/sum(N)); para(2)*(N(1)/sum(N)), para(2)*(N(2)/sum(N))];

% heterogenous mixing, varying contact rates between groups
% beta = [para(1) * ( N(1) / (N(1)+para(3)*N(2)) ), para(1) * ( para(3)*N(2) / (N(1)+para(3)*N(2)) ); para(2) * ( para(3)*N(1) / (para(3)*N(1)+N(2)) ), para(2) * ( N(2) / (para(3)*N(1)+N(2)) )];

% heterogenous mixing, varying contact rates between groups, varying pref
beta = [para(1) * ( N(1) / (N(1)+para(3)*N(2)) ), para(1) * ( para(3)*N(2) / (N(1)+para(3)*N(2)) ); para(2) * ( para(4)*N(1) / (para(4)*N(1)+N(2)) ), para(2) * ( N(2) / (para(4)*N(1)+N(2)) )];




% solve ODE given beta (along with other parameters) at initial value y1
[t,y]=ode45(@(t,y) ode_fun_two_groups(t, y, beta, T, D, N), tspan, y1);

% construct matrix with predicted case counts for each group
matrix = zeros(length(tspan),2);

for i = tspan
    matrix(i,1) = [N(1) - y(i,1)]*T(1); % reported infections in Group 1
    matrix(i,2) = [N(2) - y(i,5)]*T(2); % reported infections in Group 2
end

end