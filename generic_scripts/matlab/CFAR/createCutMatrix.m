function [full_cutidx] = createCutMatrix(detector,range_bins,doppler_bins)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% create Full Cells Under Test (CUT) matrix
        
            Ngc = detector.GuardBandSize(2);
            Ngr = detector.GuardBandSize(1);
            Ntc = detector.TrainingBandSize(2);
            Ntr = detector.TrainingBandSize(1);
            colstart = 1;
            rowstart = 1;
            
            full_cutidx = [];
      
            colend = range_bins;
            
            rowend = doppler_bins;
            
            for m = colstart:colend
                for n = rowstart:rowend
                    full_cutidx = [full_cutidx,[n;m]];
                end
            end
            
            ncutcells = size(full_cutidx,2);
            
            % show image of CUT
            cutimage = zeros(doppler_bins,range_bins);
            for k = 1:ncutcells
                cutimage(full_cutidx(1,k),full_cutidx(2,k)) = 1;
            end
            figure
            imagesc(cutimage)
            axis equal
end

