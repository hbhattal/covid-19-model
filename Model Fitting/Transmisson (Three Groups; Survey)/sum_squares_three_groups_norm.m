function sum = sum_squares_three_groups_norm(predicted, observedI, tspan)

sum = 0;

% normalized residual sum-squares calculation for reported cases
for i=tspan
    for j = 1:3
        sum = sum + 1*((predicted(i,j) - observedI(i,j))./observedI(i,j)).^2;
    end
end

end
