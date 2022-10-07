function CLEANed_cpi_matrix_array_norm = CLEAN_DSI(cross_ambg_matrix_array, self_ambg_matrix_array, ref_matrix, number_cpi, segments_per_cpi, cpi_stride, doppler_zero_padding, range_bins, max_iterations, threshold, p,range_axis,doppler_axis)
%CLEAN_DSI CLEAN Direct Signal Interference Cancellation - Post-Processing
    %   cpi_matrix_array_norm:  normalised array of range-Doppler surfaces - CAF(τ,f).
    %   ref_matrix:             segmented reference matrix 
    %   number_cpi:             number of slices per captire
    %   segments_per_cpi:       number of segments per CAF slice
    %   cpi_stride:             number of segments to stide per per CPI
    %   doppler_zero_padding:   amount of zero-padding used for Doppler FFT
    %   range_bins:             number of range bins 
    %   max_itterations:        max number of iterations operated on each slice
    %   threshold:              cutoff threshold parameter
    %   p:                      must be a positive integer less than one to avoid unwanted 
    %                           discontinuities arising from zero values in the range–Doppler surface.



%% if the input is a cell matrix
if iscell(cross_ambg_matrix_array)    
    
    %% First check there is a zero-Doppler bin - otherwise throw error
        doppler_bins = size(cross_ambg_matrix_array{1},1);
        zero_Doppler_cut_index = (doppler_bins+1)/2;
        
        if floor(zero_Doppler_cut_index)~= ceil(zero_Doppler_cut_index)
            "CAF Slices must have a zero-Doppler Doppler bin"
            return
        end
    
    
    % Initialise matrix for CLEANed data
    CLEANed_cpi_matrix_array_norm = createArrays(number_cpi, [range_bins doppler_bins]);
    
    
    %% loop through all caf slices 
        tic
        for i=1:number_cpi
         
         CLEAN_progress = 100*i/number_cpi;
         "Progress : " + ceil(CLEAN_progress) + "%"
         
            
         % populate CLEANed CAF matrix with normalised pre-CLEANed data
         CLEANed_cpi_matrix_array_norm{i} = cross_ambg_matrix_array{i}./max(cross_ambg_matrix_array{i}(:)); 
                     
         iteration = 0;   
         while iteration <= max_iterations
            iteration = iteration + 1;
            [scale_factor, delay] = max(CLEANed_cpi_matrix_array_norm{i}(zero_Doppler_cut_index,:)); % max value and delay on zero-Doppler-cut
            zero_doppler_peak = 10*log(abs(scale_factor));
    
            [max_val, index] = max(CLEANed_cpi_matrix_array_norm{i}(:));
            cpi_peak = 10*log(abs(CLEANed_cpi_matrix_array_norm{i}(8283)));
    
            relative_pwr = cpi_peak - zero_doppler_peak;
         
%             figure
%             slice = 10*log(abs(CLEANed_cpi_matrix_array_norm{i}./max(CLEANed_cpi_matrix_array_norm{i}(:)))); 
%             fig = imagesc(range_axis,doppler_axis,slice,[-inf +inf]);
%                 title("Target to Peak DSI/Clutter Relative Power : " + relative_pwr + " dB")
%                 ylabel('Doppler (Hz)')
%                 xlabel('Range (m)')   
%                 c = colorbar;
%                 c.Label.String='Norm Power (dB)';
%                 ylim([-100 100]);
%                 xlim([-inf +inf]);
%                 colormap jet;
%                 name = "Slice_" + i + "_Iteration " + iteration;
%                 saveas(fig,name,'png');
    
            if abs(scale_factor) > threshold
                 
                shifted_self_ambg = self_ambg_matrix_array{i}./max(self_ambg_matrix_array{i}(:));
                 shifted_self_ambg = circshift(shifted_self_ambg,delay-1,2);  
                 shifted_self_ambg(:,1:delay-1) = 0;
                 
                 CLEANing_surface = p * scale_factor * shifted_self_ambg;
                 CLEANed_cpi_matrix_array_norm{i} = CLEANed_cpi_matrix_array_norm{i} - CLEANing_surface;
    
            else 
                break % break from while loop to next CPI
            end
        end                                 
                
         end
    mean_process_time_per_slice = toc/number_cpi


else 
    
    %% if the matrix is just a single range doppler slice
    
    %% First check there is a zero-Doppler bin - otherwise throw error
        doppler_bins = size(cross_ambg_matrix_array,1);
        zero_Doppler_cut_index = (doppler_bins+1)/2;
        
        if floor(zero_Doppler_cut_index)~= ceil(zero_Doppler_cut_index)
            "CAF Slices must have a zero-Doppler Doppler bin"
            return;
        end
    
    
    % Initialise matrix for CLEANed data
    CLEANed_cpi_matrix_array_norm = createArrays(number_cpi, [range_bins doppler_bins]);
    relative_pwr = nan(max_iterations,1);
    
    %% loop through all caf slices 
        tic
        for i=1:number_cpi
         
         CLEAN_progress = 100*i/number_cpi;
         "Progress : " + ceil(CLEAN_progress) + "%"
         
            
         % populate CLEANed CAF matrix with normalised pre-CLEANed data
         CLEANed_cpi_matrix_array_norm = cross_ambg_matrix_array./max(cross_ambg_matrix_array(:)); 
                     
         iteration = 0;   
         while iteration <= max_iterations
            iteration = iteration + 1;
            [scale_factor, delay] = max(CLEANed_cpi_matrix_array_norm(zero_Doppler_cut_index,:)) % max value and delay on zero-Doppler-cut
            zero_doppler_peak = 10*log(abs(scale_factor)) % normalised power of zero-Doppler peak
    
            [max_val, index] = max(CLEANed_cpi_matrix_array_norm(:)) % range-Doppler surface peak
            cpi_peak = 10*log(abs(CLEANed_cpi_matrix_array_norm(8283)))

            relative_pwr(iteration) = cpi_peak - zero_doppler_peak;
          
%             figure
            slice = 10*log(abs(CLEANed_cpi_matrix_array_norm./max(CLEANed_cpi_matrix_array_norm(:)))); 
            fig = imagesc(range_axis,doppler_axis,slice,[-inf +inf]);
%                 title("Target to Peak DSI/Clutter Relative Power : " + relative_pwr(iteration) + " dB")
                ylabel('Doppler (Hz)')
                xlabel('Range (m)')   
                c = colorbar;
                c.Label.String='Norm Power (dB)';
                ylim([-100 100]);
                xlim([-inf +inf]);
                colormap jet;
%                 name = "Slice_" + i + "_Iteration " + iteration;
%                 saveas(fig,name,'png');
    
            if abs(scale_factor) > threshold
                 
                shifted_self_ambg = self_ambg_matrix_array./max(self_ambg_matrix_array(:));
                shifted_self_ambg = circshift(shifted_self_ambg,delay-1,2); % circ shift is a bad idea i think; pretty sure it results in adding DSI like artifacts  
%                 shifted_self_ambg = [zeros(size(self_ambg_matrix_array,2),delay-1), shifted_self_ambg(:,delay-1)]; % circ shift is a bad idea i think; pretty sure it results in adding DSI like artifacts  
                shifted_self_ambg(:,1:delay-1) = 0;
                 
                CLEANing_surface = p * scale_factor * shifted_self_ambg;
                CLEANed_cpi_matrix_array_norm = CLEANed_cpi_matrix_array_norm - CLEANing_surface;

%                 %plot self ambiguity surface
%                 fig = imagesc(range_axis,doppler_axis,10*log10(abs(shifted_self_ambg)),[-inf +inf]);
%                 ylabel('Doppler (Hz)')
%                 xlabel('Range (m)')   
%                 c = colorbar;
%                 c.Label.String='Norm Power (dB)';
%                 ylim([-inf inf]);
%                 xlim([-inf +inf]);
%                 colormap jet;   



            else 
                break % break from while loop to next CPI
            end
        end                                 
                
         end
    mean_process_time_per_slice = toc/number_cpi
figure
plot(relative_pwr)
xlabel('CLEAN Iteration')
ylabel('Relative Power (dB)')
grid on
grid minor 

end   

