function error = error_fun(T_test, D, N, tspan, Month, dataI, start_date, observedI)

T=[T_test,T_test];

para0 = [0.5];

predicted = @(para) pred_cases_two_groups_survey60(D, N, para, tspan, Month, dataI, start_date, T);
min_fun = @(para) sum_squares_two_groups_norm(predicted(para), observedI, tspan);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0];
ub = [1];

optimum = fmincon(min_fun, para0,A,b,Aeq,beq,lb,ub);
error = min_fun(optimum);

end