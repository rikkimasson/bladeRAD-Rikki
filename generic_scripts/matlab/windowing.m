function [outputArg1] = windowing(surv_ch, window_function)
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
            window = ones(size(surv_ch,1));
        elseif window_function == "Flat top"
            window = flattopwin(size(surv_ch,1));
        elseif window_function == "Hamming"
            window = hamming(size(surv_ch,1));
        elseif window_function == "Hann"
            window = hann(size(surv_ch,1));
%         elseif window_function == "Tukey":
%             window = signal.tukey(size(surv_ch))
        elseif window_function == "Blackman"
            window = blackman(size(surv_ch,1));
        elseif window_function == "Blackman-Harris"
            window = blackmanharris(size(surv_ch,1));
        else
            print("Window function is not identified");
        end
        outputArg1 = surv_ch .* window;

        end