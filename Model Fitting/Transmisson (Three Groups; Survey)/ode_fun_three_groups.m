function dydt = ode_fun_three_groups(t, y, beta, T, D, N)

% define compartments 
S1 = y(1);
I1R = y(2);
I1UR = y(3);
R1 = y(4);

S2 = y(5);
I2R = y(6);
I2UR = y(7);
R2 = y(8);

S3 = y(9);
I3R = y(10);
I3UR = y(11);
R3 = y(12);

% define dynamical system
dS1 = - (S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))];
dI1R = [(S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] - (1/D(1))*(I1R/T(1)) ] * T(1);
dI1UR =[(S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] - (1/D(1))*(I1R/T(1)) ] * (1-T(1));
dR1 = (1/D(1))*(I1R/T(1));

dS2 = - (S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))];
dI2R = [(S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] - (1/D(2))*(I2R/T(2)) ] * T(2);
dI2UR =[(S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] - (1/D(2))*(I2R/T(2)) ] * (1-T(2));
dR2 = (1/D(2))*(I2R/T(2));

dS3 = - (S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))];
dI3R = [(S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] - (1/D(3))*(I3R/T(3)) ] * T(3);
dI3UR =[(S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] - (1/D(3))*(I3R/T(3)) ] * (1-T(3));
dR3 = (1/D(3))*(I3R/T(3));

dydt = [dS1; dI1R; dI1UR; dR1; dS2; dI2R; dI2UR; dR2; dS3; dI3R; dI3UR; dR3];
end