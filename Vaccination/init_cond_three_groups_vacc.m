function init = init_cond_three_groups_vacc(dataI, start_date, offset, N, T)

% define relevant dates to calculate number of active cases
start_date = datetime(start_date);
start_date_minus_offset = start_date - offset;

% index rows corresponding to start_date and start_date_minus_offset from dataset
index_start_date = (dataI.Dates == start_date);
index_start_date_minus_offset = (dataI.Dates == start_date_minus_offset);

cumul_cases_start_date = dataI{index_start_date,{'Group1','Group2','Group3'}};
cumul_cases_start_date_minus_offset = dataI{index_start_date_minus_offset,{'Group1','Group2','Group3'}};

% define initial condition at t=t0
S1    = N(1) - cumul_cases_start_date(1)/T(1);
I1R   = cumul_cases_start_date(1) - cumul_cases_start_date_minus_offset(1);
I1UR  = [cumul_cases_start_date(1) - cumul_cases_start_date_minus_offset(1)]*(1-T(1))/T(1);
VS1   = 0;
VI1R  = 0; 
VI1UR = 0; 
H1    = 0;
R1    = cumul_cases_start_date_minus_offset(1);

S2    = N(2) - cumul_cases_start_date(2)/T(2);
I2R   = cumul_cases_start_date(2) - cumul_cases_start_date_minus_offset(2);
I2UR  = [cumul_cases_start_date(2) - cumul_cases_start_date_minus_offset(2)]*(1-T(2))/T(2);
VS2   = 0;
VI2R  = 0; 
VI2UR = 0; 
H2    = 0;
R2    = cumul_cases_start_date_minus_offset(2);

S3    = N(3) - cumul_cases_start_date(3)/T(3);
I3R   = cumul_cases_start_date(3) - cumul_cases_start_date_minus_offset(3);
I3UR  = [cumul_cases_start_date(3) - cumul_cases_start_date_minus_offset(3)]*(1-T(3))/T(3);
VS3   = 0;
VI3R  = 0; 
VI3UR = 0; 
H3    = 0;
R3    = cumul_cases_start_date_minus_offset(3);

init = [S1, I1R, I1UR, VS1, VI1R, VI1UR, H1, R1, S2, I2R, I2UR, VS2, VI2R, VI2UR, H2, R2, S3, I3R, I3UR, VS3, VI3R, VI3UR, H3, R3];
end