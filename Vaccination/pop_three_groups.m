function pop = pop_three_groups(age_cutoff_1, age_cutoff_2)

% Population of BC
N_T = 5.071e6; 
% Proportion of individuals in each age groups (0-9, 10-19, ..., 100+)
age_population = [0.0984	0.1060	0.1271	0.1307	0.1328	0.1526	0.1316	0.0747	0.0372	0.0087	0.0003];

% Estimate population in each group from total BC population
index1 = floor(age_cutoff_1/10);
index2 = floor(age_cutoff_2/10);

pop1 = sum(age_population(1:index1))*N_T;
pop2 = sum(age_population((index1+1):index2))*N_T;
pop3 = sum(age_population((index2+1):11))*N_T;

pop = [pop1, pop2, pop3];
end