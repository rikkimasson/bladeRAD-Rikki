function [outputArg1] = windowing(data, window_function)
%Applies the specified window function on the surveillance channel.
%         Description:
%         ------------
%             Applies the specified window function on the surveillance channel. Using the appropriate window function
%             the energy leakage of the reflected signal component in the Doppler domain can be avoided.
%             To select the desired window function chose one from the followings:
%                 - "Rectangular"
%                 - "Flat top"
%                 - "Hamming"
%                 - "Hann"
%                 - "Tukey"
%                 - "Blackman"
%                 - "Blackman-Harris"
%         Parameters:
%         -----------
%             :param: surv_ch: Surveillance channel
%             :param: window_function: (string) Name of the applied window function
%         Return values:
%         --------------
%         
        if window_function == "Rectangular"
            w = ones(size(data,1));
        elseif window_function == "Flat top"
            w = flattopwin(size(data,1));
        elseif window_function == "Hamming"
            w = hamming(size(data,1));
        elseif window_function == "Hann"
            w = hann(size(data,1));
        elseif window_function == "chebwin"
            w = chebwin(size(data,1));
        elseif window_function == "Blackman"
            w = blackman(size(data,1));
        elseif window_function == "Blackman-Harris"
            w = blackmanharris(size(data,1));
        else
            "Window function is not identified";
        end
        outputArg1 = data .* w;

        end