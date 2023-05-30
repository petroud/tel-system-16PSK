% Function to detect 16-PSK symbols using the nearest neighbor rule
% and to compute the estimated bit sequence using the inverse Gray code
function [est_X, est_bit_seq] = detect_PSK_16(Y)
    
    % Get the number of symbols (rows in Y)
    N = size(Y,1);  
    
    % Check that each symbol in Y is 2-dimensional
    if size(Y,2) ~= 2
        error('Each row of Y must be 2-dimensional for 16-PSK');
    end
    
    % Initialize the estimated symbol and bit sequence matrices
    est_X = zeros(N,2);
    est_bit_seq = zeros(N,4);
    
    % For each symbol in Y...
    for n = 1:N
        % Initialize the minimum distance to infinity
        min_dist = inf;
        
        % For each possible 16-PSK symbol...
        for m = 0:15
            % Compute the 16-PSK symbol for m
            X_m = [cos(2*pi*m/16), sin(2*pi*m/16)];
            
            % Compute the Euclidean distance from Y(n,:) to X_m
            dist = norm(Y(n,:) - X_m);
            
            % If this distance is less than the current minimum distance...
            if dist < min_dist
                % ... then update the minimum distance ...
                min_dist = dist;
                % ... and update the estimated symbol and bit sequence.
                est_X(n,:) = X_m;
                est_bit_seq(n,:) = de2bi(m, 4, 'left-msb');  % Convert integer m to 4-bit binary using Gray encoding
            end
        end
    end
end


