function val = get_multilinear_form(F_handle, varargin)
% GET_MULTILINEAR_FORM compute multilinear function with arbitrary order (non-intrusive)
%
% usage:
%   val = get_multilinear_form(F_handle, x1, x2, ..., xk)
%
% input:
%   F_handle - function handle，only one input: F_k(x)。
%   x1...xk  - input arguments, must have k arguments，where k 
%              corresponding to the order of F_handle 
%
% output:
%   val      - multilinear form: the value of F_k(x1, x2, ..., xk) 
%
% Principle:
%   base on polarization identity，and using the even and odd properties of function：
%   F(x1...xk) = 1/(2^(k-1) * k!) * sum( prod(eps) * F( sum(eps_i * x_i) ) )
%

    k = numel(varargin);
    if k < 2
        error('k should be \ge 2');
    end

    coeff = 1 / (2^(k-1) * factorial(k));

    arg_sum = varargin{1};
    for j = 2:k
        arg_sum = arg_sum + varargin{j};
    end

    total_sum = F_handle(arg_sum); 

    % Traverse the remaining 2^(k-1) -1 symbol combinations
    num_combos = 2^(k-1);
    
    for i = 1 : num_combos - 1
        % Convert i to binary, corresponding to the symbols eps_2 to eps_k.
        % 0 represents +1, 1 represents -1
        % bitget(i, j) retrieves the j-th bit.
        
        % Initialize the signed product of the current combination (eps_1 is always incremented by 1).
        sign_prod = 1;
        
        current_arg = varargin{1};
        
        for bit_idx = 1 : (k-1)
            is_negative = bitget(i, bit_idx); 
            
            input_idx = bit_idx + 1; 
            
            if is_negative
                sign_prod = -sign_prod;
                current_arg = current_arg - varargin{input_idx};
            else
                current_arg = current_arg + varargin{input_idx};
            end
        end
        
        F_val = F_handle(current_arg);
        
        if sign_prod == 1
            total_sum = total_sum + F_val;
        else
            total_sum = total_sum - F_val;
        end
    end

    val = coeff * total_sum;
end
