function num_of_symbol_errors = symbol_errors(est_X, X)
    % Compare the first column of est_X and X element-wise, and generate a logical array.
    column1_errors = est_X(:, 1) ~= X(:, 1);
    
    % Compare the second column of est_X and X element-wise, and generate a logical array.
    column2_errors = est_X(:, 2) ~= X(:, 2);
    
    % Combine the logical arrays using the logical OR operator.
    % Elements with true values indicate symbol errors.
    total_errors = column1_errors | column2_errors;
    
    % Compute the sum of true values to get the number of symbol errors.
    num_of_symbol_errors = sum(total_errors);
    
    % Return the total number of symbol errors.
end