function matrix = pred_cases_three_groups_survey_40_70(T, D, N, para, tspan, y1, Month)

% Specify which survey month (cutoff=40,70) to use as constraint
if strcmp(Month,'May')
    beta = para(1).*[1.6634746, 0.8856780, 0.1300000; 0.9661808, 1.1009038, 0.2548397; 0.6400000, 0.9200000, 0.6900000];
elseif strcmp(Month,'July')
    beta = para(1).*[1.3833813, 0.9829496, 0.0700000; 0.9155090, 1.0370659, 0.1691018; 0.2200000, 0.6400000, 0.7300000];
elseif strcmp(Month,'Sept')
    beta = para(1).*[3.2528235, 1.7652941, 0.1065882; 1.8256000, 1.7781143, 0.2797714; 0.3600000, 1.0400000, 0.5500000];
elseif strcmp(Month,'Dec')
    beta = para(1).*[1.8493976, 1.3598795, 0.2846988; 1.2988957, 1.6465644, 0.2166258; 0.6800000, 0.7400000, 0.3600000];
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