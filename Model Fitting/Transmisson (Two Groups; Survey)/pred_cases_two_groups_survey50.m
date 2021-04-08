function matrix = pred_cases_two_groups_survey50(T, D, N, para, tspan, y1, Month)

% Specify which survey month (cutoff=50) to use as constraint
if strcmp(Month,'May')
    beta = para(1)*[1.938461538, 0.711481481; 1.012470588, 1.195088235];
elseif strcmp(Month,'July')
    beta = para(1)*[1.866871795, 0.607846154; 0.848194444, 0.977291667];
elseif strcmp(Month,'Sept')
    beta = para(1)*[3.786944444, 1.286666667; 1.542674419, 1.328313953];
elseif strcmp(Month,'Dec')
    beta = para(1).*[2.202071429, 1.196642857; 1.532119205, 1.181125828];
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