function [dets, thres, noise_pwr] = cfarSlicesDetections(range_doppler_slices,detector,cutidx,range_lookup,Doppler_lookup)
%CFARSLICES loop through matrix and operate CFAR detector
    number_cpi = size(range_doppler_slices,2); 
    dets = createArrays(number_cpi,[1 1]);         
    thres = createArrays(number_cpi,[1 1]); 
    noise_pwr = createArrays(number_cpi,[1 1]); 

    for i=1:number_cpi
        i
        % CFAR the CAF to get detections array dets
        [detections,thres{i},noise_pwr{i}] = detector(abs(range_doppler_slices{i}),cutidx);
        % convert detection cell indexes to range and Doppler values
        range = range_lookup(detections(2,:));
        doppler = Doppler_lookup(detections(1,:));
        dets{i} = [ range; doppler ];
    end

end

