function prbs_seq = generate_prbs(n, init_state)
% GENERATE_PRBS Generate PRBS sequence based on polynomial X^11 + X^2 + 1
%
% Inputs:
%   n          - Number of bits to generate
%   init_state - (Optional) Initial 11-bit state as vector [default: all ones]
%
% Output:
%   prbs_seq   - PRBS sequence of length n (column vector)
%
% Example:
%   seq = generate_prbs(100);              % Generate 100 bits
%   seq = generate_prbs(50, [1 0 1 0 1 0 1 0 1 0 1]); % Custom init state

    % Default initial state (all ones) if not provided
    if nargin < 2
        init_state = ones(1, 11);
    end
    
    % Ensure init_state is the correct length
    if length(init_state) ~= 11
        error('Initial state must be 11 bits long');
    end
    
    % Initialize shift register (convert to column vector)
    shift_reg = init_state(:);
    
    % Preallocate output
    prbs_seq = zeros(n, 1);
    
    % Generate PRBS sequence
    % Polynomial: X^11 + X^2 + 1
    % Feedback taps at positions 11 and 2
    for i = 1:n
        % Output is the last bit in the shift register
        prbs_seq(i) = shift_reg(end);
        
        % Calculate feedback: XOR of tap positions 11 and 2
        feedback = xor(shift_reg(11), shift_reg(9));
        
        % Shift register right and insert feedback at the beginning
        shift_reg = [feedback; shift_reg(1:end-1)];
    end
    
end