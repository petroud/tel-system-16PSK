function num_of_bit_errors = bit_errors(est_bit_seq, b)
    % Check the dimensions of the input bit sequences
    if size(est_bit_seq) ~= size(b)
        error('Dimensions of the estimated bit sequence and the true bit sequence must match');
    end
    
    % Compute the bit errors
    bit_errors = est_bit_seq ~= b;
    
    % Sum all the bit errors
    num_of_bit_errors = sum(bit_errors(:));
end