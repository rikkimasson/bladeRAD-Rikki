function [lower, upper] = findClosestNeighbors(array, n)
    % findClosestNeighbors - Find closest values on either side of n
    %
    % Inputs:
    %   array - sorted array of numbers (ascending order)
    %   n     - target number
    %
    % Outputs:
    %   lower - closest value less than or equal to n
    %   upper - closest value greater than or equal to n
    
    % Validate inputs
    if ~issorted(array)
        error('Array must be sorted in ascending order');
    end
    
    if n < array(1) || n > array(end)
        error('n must be between the smallest and largest values in the array');
    end
    
    % Find the index where n would be inserted to maintain sorted order
    idx = find(array >= n, 1, 'first');
    
    % Handle edge cases
    if isempty(idx)
        % n is larger than all elements (shouldn't happen due to validation)
        lower = array(end);
        upper = array(end);
    elseif array(idx) == n
        % Exact match found
        lower = array(idx);
        upper = array(idx);
    elseif idx == 1
        % n is smaller than first element (shouldn't happen due to validation)
        lower = array(1);
        upper = array(1);
    else
        % Normal case: n is between two elements
        lower = (idx - 1);
        upper = (idx);
    end
end
