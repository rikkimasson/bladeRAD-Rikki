function [dets, thres, noise_pwr] = cfarSlicesDetections(range_doppler_slices,detector,cutidx,range_lookup,Doppler_lookup)
%CFARSLICES loop through matrix and operate CFAR detector
    number_cpi = size(range_doppler_slices,2); 
    dets = createArrays(number_cpi,[1 1]);         
    thres = createArrays(number_cpi,[1 1]); 
    noise_pwr = createArrays(number_cpi,[1 1]);
    numCUT = size(range_doppler_slices{1},1)*size(range_doppler_slices{1},2);
    trninds = createArrays(numCUT, [1 1]);

    for i=1:number_cpi
        tic
        % CFAR the CAF to get detections array dets
        [detections,thres{i},noise_pwr{i}] = detector(abs(range_doppler_slices{i}),cutidx);
        % convert detection cell indexes to range and Doppler values
        range = range_lookup(detections(2,:));
        doppler = Doppler_lookup(detections(1,:));
        dets{i} = [ range; doppler ];
        duration = toc;
        time_remaining = toc*(number_cpi-i);
        time_remaining = "Time remaining : " + time_remaining;
        cpis = "CPI : " + i + " of " + number_cpi;
        display(cpis)
        display(time_remaining)
    end

end

