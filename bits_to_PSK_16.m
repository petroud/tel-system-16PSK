function X = bits_to_PSK_16(bit_seq)
    N = size(bit_seq,1);  % number of rows
    if size(bit_seq,2) ~= 4
        error('Each row of bit sequence must be 4 bits long for 16-PSK');
    end
    X = zeros(N,2);
    for n = 1:N
        bits = bit_seq(n,:);
        m = bi2de(bits,'left-msb');  % Convert bits to integer using Gray encoding
        X(n,:) = [cos(2*pi*m/16), sin(2*pi*m/16)];
    end
end