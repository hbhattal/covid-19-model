function sum = sum_squares_two_groups(predictedI, observedI, tspan)

sum = 0;

% residual sum-squares calculation for reported cases
for i=tspan
    for j = 1:2
        sum = sum + 1*(predictedI(i,j) - observedI(i,j)).^2;
    end
end

end
