function matrix = pred_cases_three_groups_survey_50_70(T, D, N, para, tspan, y1, Month)

% Specify which survey month (cutoff=50,70) to use as constraint
if strcmp(Month,'May')
    beta = para(1).*[1.9384615, 0.5487179, 0.1627635; 1.0578947, 0.8614474, 0.2673684; 0.9200000, 0.6400000, 0.6900000];
elseif strcmp(Month,'July')
    beta = para(1).*[1.8668718, 0.5234872, 0.0843590; 0.9903604, 0.7112613, 0.1938739; 0.3700000, 0.4900000, 0.7300000];
elseif strcmp(Month,'Sept')
    beta = para(1).*[3.7869444, 1.0967361, 0.1899306; 1.9301724, 1.1210345, 0.2643966; 0.7400000, 0.6600000, 0.5500000];
elseif strcmp(Month,'Dec')
    beta = para(1).*[2.2020714, 0.9505000, 0.2461429; 1.8047170, 1.0737736, 0.2309434; 0.8900000, 0.5300000, 0.3600000];
end


% solve ODE given beta (along with other parameters) at initial value y1
[t,y]=ode45(@(t,y) ode_fun_three_groups(t, y, beta, T, D, N), tspan, y1);

% construct matrix with predicted case counts for each group
matrix = zeros(length(tspan),3);

for i = tspan
    matrix(i,1) = [N(1) - y(i,1)]*T(1); % reported infections in Group 1
    matrix(i,2) = [N(2) - y(i,5)]*T(2); % reported infections in Group 2
    matrix(i,3) = [N(3) - y(i,9)]*T(3); % reported infections in Group 2
end

end