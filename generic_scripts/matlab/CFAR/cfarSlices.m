function [cfar_range_doppler_slices] = cfarSlices(range_doppler_slices,detector,cutidx)
%CFARSLICES loop through matrix and operate CFAR detector

  number_cpi = size(range_doppler_slices,2); 
  range_bins = size(range_doppler_slices{1},1);
  doppler_bins = size(range_doppler_slices{1},2);
  % CAF slices after CFAR (binary detections)
    cfar_range_doppler_slices = createArrays(number_cpi, [range_bins doppler_bins]);         

    for i=1:number_cpi
        i
        % CFAR the CAF to get detections array dets
        [dets,th1] = detector(abs(range_doppler_slices{i}),cutidx);
        % reshape array to matrix
        caf = reshape(dets,[range_bins],[doppler_bins]);
        cfar_range_doppler_slices{i} = caf;

    end

end

