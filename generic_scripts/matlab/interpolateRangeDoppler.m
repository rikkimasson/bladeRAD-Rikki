function [interRangeDoppler] = interpolateRangeDoppler(rangeDopplerSlices,rangeFactor,dopplerFactor)
%INTERPOLATERANGEDOPPLER Summary of this function goes here
%   Detailed explanation goes here

%% if the input is a cell matrix
if iscell(rangeDopplerSlices)    
    
    nSurfs = size(rangeDopplerSlices,2);
    [nDopBins nRangeBins] = size(rangeDopplerSlices{1});

    interpnDopBins = round(nDopBins * dopplerFactor);
    interpnRangeBins = round(nRangeBins * rangeFactor);

     if mod(interpnDopBins,2) == 0
      interpnDopBins = interpnDopBins + 1;
     end


    [x, y] = meshgrid(1:nRangeBins, 1:nDopBins);
    [xq, yq] = meshgrid(linspace(1, nRangeBins, interpnRangeBins), linspace(1, nDopBins, interpnDopBins));



    % Initialise matrix for CLEANed data
    interRangeDoppler = createArrays(nSurfs, [interpnDopBins interpnRangeBins]);
    
    
    %% loop through all range Doppler Surfaces        
        for i=1:nSurfs
            interRangeDoppler{i} = interp2(x, y, rangeDopplerSlices{i}, xq, yq,'makima');
        end


end

