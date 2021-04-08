function matrix = pred_cases_two_groups_survey70(T, D, N, para, tspan, y1, Month)

% Specify which survey month (cutoff=70) to use as constraint
if strcmp(Month,'May')
    beta = para(1).*[2.26357513, 0.203955095; 1.56, 0.69];
elseif strcmp(Month,'July')
    beta = para(1).*[2.140522876, 0.124084967; 0.86, 0.73];
elseif strcmp(Month,'Sept')
    beta = para(1).*[4.066115385, 0.223153846; 1.4, 0.55];
elseif strcmp(Month,'Dec')
    beta = para(1).*[3.034471545, 0.239593496; 1.42, 0.36];
end


% solve ODE given beta (along with other parameters) at initial value y1
[t,y]=ode45(@(t,y) ode_fun_two_groups(t, y, beta, T, D, N), tspan, y1);

% construct matrix with predicted case counts for each group
matrix = zeros(length(tspan),2);

for i = tspan
    matrix(i,1) = [N(1) - y(i,1)]*T(1); % reported infections in Group 1
    matrix(i,2) = [N(2) - y(i,5)]*T(2); % reported infections in Group 2
end

end