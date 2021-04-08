function dydt = ode_fun_two_groups(t, y, beta, T, D, N)

% define compartments 
S1 = y(1);
I1R = y(2);
I1UR = y(3);
R1 = y(4);

S2 = y(5);
I2R = y(6);
I2UR = y(7);
R2 = y(8);

% define dynamical system
dS1 = - (S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2))];
dI1R = [(S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2))] - (1/D(1))*(I1R/T(1)) ] * T(1);
dI1UR =[(S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2))] - (1/D(1))*(I1R/T(1)) ] * (1-T(1));
dR1 = (1/D(1))*(I1R/T(1));

dS2 = - (S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2))];
dI2R = [(S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2))] - (1/D(2))*(I2R/T(2)) ] * T(2);
dI2UR =[(S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2))] - (1/D(2))*(I2R/T(2)) ] * (1-T(2));
dR2 = (1/D(2))*(I2R/T(2));

dydt = [dS1; dI1R; dI1UR; dR1; dS2; dI2R; dI2UR; dR2];
end