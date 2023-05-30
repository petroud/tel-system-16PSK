function num_of_symbol_errors = symbol_errors(est_X, X)
    % Compute symbol errors
    num_of_symbol_errors = sum(est_X(:,1) ~= X(:,1) | est_X(:,2) ~= X(:,2));
end