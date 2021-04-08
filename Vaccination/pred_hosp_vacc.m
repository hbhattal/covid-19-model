function y = pred_hosp_vacc(para, T, D, N, V_T, V_PI, V_PS, HF, days_vaccinating, vacc_per_day, vacc_prop, tspan, y1, Month, strategy)

% Contact matrices for ages: 0-49, 50-69, 70+
if strcmp(Month,'May')
    beta = para(1).*[1.9384615, 0.5487179, 0.1627635; 1.0578947, 0.8614474, 0.2673684; 0.9200000, 0.6400000, 0.6900000];
elseif strcmp(Month,'July')
    beta = para(1).*[1.8668718, 0.5234872, 0.0843590; 0.9903604, 0.7112613, 0.1938739; 0.3700000, 0.4900000, 0.7300000];
elseif strcmp(Month,'Sept')
    beta = para(1).*[3.7869444, 1.0967361, 0.1899306; 1.9301724, 1.1210345, 0.2643966; 0.7400000, 0.6600000, 0.5500000];
elseif strcmp(Month,'Dec')
    beta = para(1).*[2.2020714, 0.9505000, 0.2461429; 1.8047170, 1.0737736, 0.2309434; 0.8900000, 0.5300000, 0.3600000];
end

% Run ODE function 
[t,y]=ode45(@(t,y) ode_fun_three_groups_vacc(t, y, beta, T, D, N, V_T, V_PI, V_PS, HF, days_vaccinating, vacc_per_day, vacc_prop, strategy), tspan, y1);

y_diff = y(length(tspan), :) - y1; 
new_hosp = [y_diff(7),  y_diff(15), y_diff(23)];
total_new_hosp = sum(new_hosp);

end