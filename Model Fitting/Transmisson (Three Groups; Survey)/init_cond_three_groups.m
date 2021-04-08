function init = init_cond_three_groups(dataI, start_date, offset, N, T)

% define relevant dates to calculate number of active cases
start_date = datetime(start_date);
start_date_minus_offset = start_date - offset;

% index rows corresponding to start_date and start_date_minus_offset from dataset
index_start_date = (dataI.Dates == start_date);
index_start_date_minus_offset = (dataI.Dates == start_date_minus_offset);

cumul_cases_start_date = dataI{index_start_date,{'Group1','Group2','Group3'}};
cumul_cases_start_date_minus_offset = dataI{index_start_date_minus_offset,{'Group1','Group2','Group3'}};

% define initial condition <S1, I1R, I1UR, R1, S2, I2R, I2UR, R2, S3, I3R, I3UR, R3> at t=t0
S1 = N(1) - cumul_cases_start_date(1)/T(1);
I1R = cumul_cases_start_date(1) - cumul_cases_start_date_minus_offset(1);
I1UR = [cumul_cases_start_date(1) - cumul_cases_start_date_minus_offset(1)]*(1-T(1))/T(1);
R1 = cumul_cases_start_date_minus_offset(1);

S2 = N(2) - cumul_cases_start_date(2)/T(2);
I2R = cumul_cases_start_date(2) - cumul_cases_start_date_minus_offset(2);
I2UR = [cumul_cases_start_date(2) - cumul_cases_start_date_minus_offset(2)]*(1-T(2))/T(2);
R2 = cumul_cases_start_date_minus_offset(2);

S3 = N(3) - cumul_cases_start_date(3)/T(3);
I3R = cumul_cases_start_date(3) - cumul_cases_start_date_minus_offset(3);
I3UR = [cumul_cases_start_date(3) - cumul_cases_start_date_minus_offset(3)]*(1-T(3))/T(3);
R3 = cumul_cases_start_date_minus_offset(3);

init = [S1, I1R, I1UR, R1, S2, I2R, I2UR, R2, S3, I3R, I3UR, R3];

end