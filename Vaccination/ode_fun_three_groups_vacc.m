function  dydt = ode_fun_three_groups_vacc(t, y, beta, T, D, N, V_T, V_PI, V_PS, HF, days_vaccinating, vacc_per_day, vacc_prop, strategy)

% Define compartments 
S1    = y(1);
I1R   = y(2);
I1UR  = y(3);
VS1   = y(4);
VI1R  = y(5);
VI1UR = y(6);
H1    = y(7);
R1    = y(8);

S2    = y(9);
I2R   = y(10);
I2UR  = y(11);
VS2   = y(12);
VI2R  = y(13);
VI2UR = y(14);
H2    = y(15);
R2    = y(16);

S3    = y(17);
I3R   = y(18);
I3UR  = y(19);
VS3   = y(20);
VI3R  = y(21);
VI3UR = y(22);
H3    = y(23);
R3    = y(24);

% Ensure susceptible population never negative
if S1 <= 0
    S1 = 0;
end
if S2<=0
    S2=0;
end
if S3<=0
    S3=0;
end

% Ensure vaccinated infective population never negative
if VS1<=0
    VS1=0;
end
if VS2<=0
    VS2=0;
end
if VS3<=0
    VS3=0;
end

% Designate strategy
if string(strategy) == 'oldest'
    if S3>0
        vacc_prop = [0,0,1];
    elseif S2>0
        vacc_prop = [0,1,0];
    elseif S1>0
        vacc_prop = [1,0,0];
    end
elseif string(strategy) == 'youngest'
    if S1>0
        vacc_prop = [1,0,0];
    elseif S2>0
        vacc_prop = [0,1,0];
    elseif S3>0
        vacc_prop = [0,0,1];
    end
end



% Ensure not vaccinating population with no susceptible
if S1<=0 
    vacc_prop(1) = 0;
end
if S2<=0 
    vacc_prop(2) = 0;
end
if S3<=0 
    vacc_prop(3) = 0;
end

% Time-dependent vaccination 
if t <= (days_vaccinating + 1)
    V = [vacc_per_day.*vacc_prop(1), vacc_per_day.*vacc_prop(2), vacc_per_day.*vacc_prop(3)];
else
    V = [0, 0, 0];
end


% Define dynamical system
dS1    = - [ (S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (S1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] ] - V(1);
dI1R   =   [ (S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (S1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] - (1/D(1))*(I1R/T(1)) ] * T(1);
dI1UR  =   [ (S1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (S1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] - (1/D(1))*(I1R/T(1)) ] * (1-T(1));
dVS1   = - (1-V_PI) * [ (VS1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (VS1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] ] + V(1);
dVI1R  =   (1-V_PI) * [ (VS1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (VS1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] - (1/D(1))*(VI1R/T(1)) ] * T(1);
dVI1UR =   (1-V_PI) * [ (VS1/N(1)) * [beta(1,1)*(I1R/T(1)) + beta(2,1)*(I2R/T(2)) + beta(3,1)*(I3R/T(3))] + (VS1/N(1)) * (1-V_T) * [beta(1,1)*(VI1R/T(1)) + beta(2,1)*(VI2R/T(2)) + beta(3,1)*(VI3R/T(3))] - (1/D(1))*(VI1R/T(1)) ] * (1-T(1));
dH1    = HF(1)/D(1) * ( I1R + (1-V_PS)*VI1R );
dR1    = 1/D(1) * (I1R/T(1) + VI1R/T(1)) - HF(1)/D(1) * (I1R + (1-V_PS)*VI1R);

dS2    = - [ (S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (S2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] ] - V(2);
dI2R   =   [ (S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (S2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] - (1/D(2))*(I2R/T(2)) ] * T(2);
dI2UR  =   [ (S2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (S2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] - (1/D(2))*(I2R/T(2)) ] * (1-T(2));
dVS2   = - (1-V_PI) * [ (VS2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (VS2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] ] + V(2);
dVI2R  =   (1-V_PI) * [ (VS2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (VS2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] - (1/D(2))*(VI2R/T(2)) ] * T(2);
dVI2UR =   (1-V_PI) * [ (VS2/N(2)) * [beta(1,2)*(I1R/T(1)) + beta(2,2)*(I2R/T(2)) + beta(3,2)*(I3R/T(3))] + (VS2/N(2)) * (1-V_T) * [beta(1,2)*(VI1R/T(1)) + beta(2,2)*(VI2R/T(2)) + beta(3,2)*(VI3R/T(3))] - (1/D(2))*(VI2R/T(2)) ] * (1-T(2));
dH2    = HF(2)/D(2) * ( I2R + (1-V_PS)*VI2R );
dR2    = 1/D(2) * (I2R/T(2) + VI2R/T(2)) - HF(2)/D(2) * (I2R + (1-V_PS)*VI2R);

dS3    = - [ (S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (S3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] ] - V(3);
dI3R   =   [ (S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (S3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] - (1/D(3))*(I3R/T(3)) ] * T(3);
dI3UR  =   [ (S3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (S3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] - (1/D(3))*(I3R/T(3)) ] * (1-T(3));
dVS3   = - (1-V_PI) * [ (VS3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (VS3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] ] + V(3);
dVI3R  =   (1-V_PI) * [ (VS3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (VS3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] - (1/D(3))*(VI3R/T(3)) ] * T(3);
dVI3UR =   (1-V_PI) * [ (VS3/N(3)) * [beta(1,3)*(I1R/T(1)) + beta(2,3)*(I2R/T(2)) + beta(3,3)*(I3R/T(3))] + (VS3/N(3)) * (1-V_T) * [beta(1,3)*(VI1R/T(1)) + beta(2,3)*(VI2R/T(2)) + beta(3,3)*(VI3R/T(3))] - (1/D(3))*(VI3R/T(3)) ] * (1-T(3));
dH3    = HF(3)/D(3) * ( I3R + (1-V_PS)*VI3R );
dR3    = 1/D(3) * (I3R/T(3) + VI3R/T(3)) - HF(3)/D(3) * (I3R + (1-V_PS)*VI3R);

dydt = [dS1; dI1R; dI1UR; dVS1; dVI1R; dVI1UR; dH1; dR1; dS2; dI2R; dI2UR; dVS2; dVI2R; dVI2UR; dH2; dR2; dS3; dI3R; dI3UR; dVS3; dVI3R; dVI3UR; dH3; dR3];
end