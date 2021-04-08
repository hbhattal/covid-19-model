function matrix = pred_cases_two_groups_survey60(T, D, N, para, tspan, y1, Month)

% Specify which survey month (cutoff=60) to use as constraint
if strcmp(Month,'May')
    beta = para(1).*[2.178040541, 0.386036036; 1.235465587, 0.959878543];
elseif strcmp(Month,'July')
    beta = para(1).*[1.912062257, 0.341984436; 1.067804878, 0.958414634];
elseif strcmp(Month,'Sept')
    beta = para(1).*[4.44574359, 0.506615385; 1.399338843, 0.738677686];
elseif strcmp(Month,'Dec')
    beta = para(1).*[2.768723404, 0.602234043; 1.661165049, 0.783300971];
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