

% Opens the file gpsLog.txt with read access.
    fileID = fopen('/media/piers/T7/15_09_2022_farm/gps_ground_truth _data/geode_logger/Nmea_20220915_114115.txt','r');
%     fileID = fopen('/media/piers/T7/15_09_2022_farm/gps_ground_truth _data/iphone_logger/2022-09-15 15_07_51.nmea','r');
%     fileID = fopen('/media/piers/T7/15_09_2022_farm/gps_ground_truth _data/iphone_logger/2022-09-15 16_10_52.nmea','r');
% Read the text file.
    gpsData = fscanf(fileID,'%c');
% Create NMEA parser
    pnmea = nmeaParser('MessageId',["GGA"]);
% Parse the NMEA Data.
    ggaData = pnmea(gpsData);
% Create NMEA parser
    pnmea = nmeaParser('MessageId',["RMC"]);
% Parse the NMEA Data.
    rmcData = pnmea(gpsData);


%% Drone GPS Ground Truth 

% drone_file_name = 'DJIFlightRecord_2019-12-11_(09-55-56).xls'; % Drone .csv data file name
time_shift =  - 16; %time shift in 10ths of a second
range_shift = 2; %meters

%DVB-T Tx Cordinates
Tx_lat = 52.130139;
Tx_lon = -0.241389;
Tx_height = 234 
[x y z] = latlonalt_xyz(Tx_lon,Tx_lat,Tx_height)
sandyheath(1) = x; 
sandyheath(2) = y;    
sandyheath(3) = z;


geode1.lat = vertcat(ggaData.Latitude);
geode1.lon = vertcat(ggaData.Longitude);
geode1.alt = vertcat(ggaData.Altitude);
geode1.time_utc = vertcat(ggaData.UTCTime)
geode1.SoG = vertcat(rmcData.GroundSpeed);




% Plot the position in geographic coordinates
figure
geoplot(geode1.lat, geode1.lon,'Marker',"*",'MarkerSize',3, ...
    "Color",'red','MarkerFaceColor','red');
% Selects the basemap
geobasemap 'streets';

% Plot location of passive radar Tx and Radar
figure
geoplot(geode1.lat(1),geode1.lon(1),'r-*',...
    Tx_lat, Tx_lon,'b*:');
geobasemap 'topographic';
text(Tx_lat,Tx_lon,'DVB-T Tower',...
    'VerticalAlignment','middle','HorizontalAlignment','right');

figure
geoplot([geode1.lat(1) Tx_lat],[geode1.lon(1) Tx_lon],'black-');
hold on
geoplot([geode1.lat(1) Tx_lat],[geode1.lon(1) Tx_lon],'r-*');
geobasemap 'topographic';
legend('Target Motion','Bistatic Baseline: 8.70 km')

text(Tx_lat,Tx_lon,'DVB-T Tower',...
    'VerticalAlignment','middle','HorizontalAlignment','right');

%% Convert lat long Altitude to [x y z] 
%   for converting from latitude, longitude, altitude (above reference ellipsoid) to Cartesian ECEF
xyz = zeros(size(ggaData,1),3);
measured_gps_range = nan(size(ggaData,1),1);
for i=1:size(ggaData,1)
   [x y z] = latlonalt_xyz(geode1.lon(i),geode1.lat(i),geode1.alt(i));
   xyz(i,1) = x;
   xyz(i,2) = y;
   xyz(i,3) = z;
   measured_gps_range(i) = straight_line_dist(radar_position_xyz(1),radar_position_xyz(2),radar_position_xyz(3),x,y,z);
end
geode1.target_range = measured_gps_range
% geode1.baseline = baseline
% mean_pos = [mean(xyz(:,1)),mean(xyz(:,2)),mean(xyz(:,3))]
figure
plot(xyz(:,1)-radar_position_xyz(1))
hold on 
plot(xyz(:,2)-radar_position_xyz(2))
hold on
plot(xyz(:,3)-radar_position_xyz(3))
grid on 
title('Survey Error of Radar - Cartesian ECEF')
xlabel('Time')
ylabel('Deviation from Mean (m)')
legend('x','y','z')
hold on 

figure
plot(geode1.target_range)
% title('GPS Measured Bistatic Range')
xlabel('Time (s)')
ylabel('Target Range (m)')
grid on 
hold on 

figure
plot(measured_gps_range)
title('GPS Measured Bistatic Range')
xlabel('Time')
ylabel('Bistatic Range (m)')
grid on 
hold on 

baseline = straight_line_dist(radar_position_xyz(1),radar_position_xyz(2),radar_position_xyz(3),sandyheath(1),sandyheath(2),sandyheath(3)) %baseline between N0 & Tx
baseline/1000

% 
% 
% 
% 
% 
% [x0, y0, z0] = latlonalt_xyz(N0_lon,N0_lat,N0_height); % covere_file_name];
% drone_gps_data = Import_gps_log(drone_data_file); % load gps data from .csv file format is [Time(s) Lat Long Alt(m)]
% 
% %%% Node0 coordinates to ECEF xyz
% 
% %Calculate Bistatic Baseline Distance 
% 
% baseline = straight_line_dist(x0,y0,z0,x1,y1,z1); %baseline between N0 & Tx
% drone_data_file = [drone_data_dir,date,slash,dron
% %calculate distance from drone to Node 0 
% 
% range = zeros(size(drone_gps_data,1),1);
% for i = 1:size(drone_gps_data,1)
%     
%     drone_lat = drone_gps_data(i,2);
%     drone_lon = drone_gps_data(i,3);
%     drone_alt = drone_gps_data(i,4) + N0_height; %account for drones height being not refenced to WGS84, use altutude of monostatic node to correct drone altitudes
%     [xd, yd, zd] = latlonalt_xyz(drone_lon,drone_lat,drone_alt);
%     range(i) = straight_line_dist(x0,y0,z0,xd,yd,zd);
% end
% 
% time_range = [drone_gps_data(:,1),range];
% 
% 
% % 2019_12_11_10_4_32
% 
% % Capture start time 10_04_32 
% % Capture duration 90000 pulses (90 Sec)
% % Drone GPS Log start time 09_55_56
% % Delta 8min 36s ; drone time base 1/10 second
% 
% delta_sec = (8 * 60) + 36 + time_shift; %delay between deon log and radar capture start
% delta_samp = (delta_sec * 10); %start location of Drone in radar capture
% end_samp = delta_samp + 899; %90 second capture
% drone_flight_path = time_range(delta_samp:end_samp,:);
% % plot(drone_flight_path(:,2));
% min_drone_range = min(drone_flight_path(:,2));
% max_drone_range = max(drone_flight_path(:,2));
% min_rbin = fix(min_drone_range * 2 / 3.3)-100;
% max_rbin = ceil(max_drone_range * 2 / 3.3)+100;
% 
% %interpolate to get location for each pulse
% new_sample_points = 1:PULSES;
% old_sample_points = linspace(1,PULSES,900);
% interp_drone_flight_path = interp1(old_sample_points,drone_flight_path(:,2),new_sample_points,'linear');
% interp_drone_flight_path_rbin = (interp_drone_flight_path*2/3.3)+range_shift; % covert from meters to range bins
%                
%         % Range Limited Plot 
%                 figure
%                 FS = 300e6;
%                 T = 1/FS;
%                 time_axis = linspace(0,Data_size_seconds,size(Data_matched,1));
%                 RTI_plot= 10*log10(abs(Data_matched./max(Data_matched(:))));
%                 RTI_plot_lim = RTI_plot(:,min_rbin:max_rbin);
%                 imagesc(min_rbin:max_rbin,time_axis,RTI_plot_lim,[-50,0]);               
%                 hold on 
%                 plot(interp_drone_flight_path_rbin,time_axis,'rd', 'MarkerSize', 2);
%                 grid on            
%                 colorbar
%                 ylabel('Time (Sec)')
%                 xlabel('Range Bin')            
%                 title(mytitle,'Interpreter', 'none');
%             %                 saveas(gcf,[file_name,'\',mytitle,'_RTI_RLIM.png']);
%                 pause(1)
% 
% 
% 
% 
