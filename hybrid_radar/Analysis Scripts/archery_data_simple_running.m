clear
clc
close all
tic
% 1. Respond in the command window when prompted in order to point MATLAB
% to the correct directory
% prompt = "Which computer is running the script? 1=Rikki's Mac, 2=Jacob's Laptop, 3=Rikki UCL \n"
pc = 3;%input(prompt);
addpath("C:/Users/rw/Documents/RADCOM/misc_helper_functions");

local_save_directory = "D:/passive_archery_data/run_5_32x_small/";
local_save_directory = "C:/Users/rw/Downloads/run_5_32x_small/run_5_32x_small/";
repo_directory = "C:\Users\rw\Documents\";


% 2. Processing Flags - choose whether to
process_active_flag = false;
process_passive_flag = true;
ground_truth_flag = false;

% 3. Select the experiment you wish to process


% adds paths to generic functions

addpath(repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab',...
repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab\CFAR\',...
repo_directory + '\bladeRAD-Rikki\generic_scripts',...
repo_directory + '\bladeRAD-Rikki\generic_scripts\ref_signals\')
  


%% Load .mat file containing experiment setup parameters

    % mat_file_name = local_save_directory + i + bracket + "Experimental Configuration.mat";
    % load(mat_file_name);
    exp_dir = local_save_directory ;


%% Passive Processing

    if process_passive_flag == true
        passive.min_range = 2; % maximum number of range bins - xcorr shifts.
        passive.max_range = 20; % maximum number of range bins - xcorr shifts.
        passive.range_zero_padding = 1; % 1 = none, 2 = 100%
        passive.td_corr = false; % true = time domain xcorr; false = freq domain xcorr
        passive.seg_s = 2000; % number of segments per second - analagos to PRF.
        passive.seg_percent = 100; % % of segment used for xcorr.
        
        % load signal and split ref and sur
        file_location = local_save_directory + 'passive_data_test_1/';
        load(file_location+'passive_ref_output');
        data=double(data);
        ref_channel  = data(1:2:end)+1j .*data(2:2:end);

        file_location = local_save_directory + 'passive_data_test_1/';
        load(file_location+'passive_surv_output');
        data=double(data);
        sur_channel  = data(1:2:end)+1j .*data(2:2:end);

        % write_binary_to_text('passive_ref',ref_channel(1:2000000),16,0);
        % 
        % write_binary_to_text('passive_surv',sur_channel(1:2000000),16,0);
    %%
permanant_carriers= [0 48 54 87 141 156 192 ...
201 255 279 282 333 432 450 ...
483 525 531 618 636 714 759 ...
765 780 804 873 888 918 939 ...
942 969 984 1050 1101 1107 1110 ...
1137 1140 1146 1206 1269 1323 1377 ...
1491 1683 1704 1752 1758 1791 1845 ...
1860 1896 1905 1959 1983 1986 2037 ...
2136 2154 2187 2229 2235 2322 2340 ...
2418 2463 2469 2484 2508 2577 2592 ...
2622 2643 2646 2673 2688 2754 2805 ...
2811 2814 2841 2844 2850 2910 2973 ...
3027 3081 3195 3387 3408 3456 3462 ...
3495 3549 3564 3600 3609 3663 3687 ...
3690 3741 3840 3858 3891 3933 3939 ...
4026 4044 4122 4167 4173 4188 4212 ...
4281 4296 4326 4347 4350 4377 4392 ...
4458 4509 4515 4518 4545 4548 4554 ...
4614 4677 4731 4785 4899 5091 5112 ...
5160 5166 5199 5253 5268 5304 5313 ...
5367 5391 5394 5445 5544 5562 5595 ...
5637 5643 5730 5748 5826 5871 5877 ...
5892 5916 5985 6000 6030 6051 6054 ...
6081 6096 6162 6213 6219 6222 6249 ...
6252 6258 6318 6381 6435 6489 6603 ...
6795 6816];

permanant_carriers=permanant_carriers+1;

%% mapping from 64qam to bits

map_64qam(1,1)='100000';
map_64qam(1,2)='100010';
map_64qam(1,3)='101010';
map_64qam(1,4)='101000';
map_64qam(1,5)='001000';
map_64qam(1,6)='001010';
map_64qam(1,7)='000010';
map_64qam(1,8)='000000';

map_64qam(2,1)='100001';
map_64qam(2,2)='100011';
map_64qam(2,3)='101011';
map_64qam(2,4)='101001';
map_64qam(2,5)='001001';
map_64qam(2,6)='001011';
map_64qam(2,7)='000011';
map_64qam(2,8)='000001';

map_64qam(3,1)='100101';
map_64qam(3,2)='100111';
map_64qam(3,3)='101111';
map_64qam(3,4)='101101';
map_64qam(3,5)='001101';
map_64qam(3,6)='001111';
map_64qam(3,7)='000111';
map_64qam(3,8)='000101';

map_64qam(4,1)='100100';
map_64qam(4,2)='100110';
map_64qam(4,3)='101110';
map_64qam(4,4)='101100';
map_64qam(4,5)='001100';
map_64qam(4,6)='001110';
map_64qam(4,7)='000110';
map_64qam(4,8)='000100';

map_64qam(5,1)='110100';
map_64qam(5,2)='110110';
map_64qam(5,3)='111110';
map_64qam(5,4)='111100';
map_64qam(5,5)='011100';
map_64qam(5,6)='011110';
map_64qam(5,7)='010110';
map_64qam(5,8)='010100';

map_64qam(6,1)='110101';
map_64qam(6,2)='110111';
map_64qam(6,3)='111111';
map_64qam(6,4)='111101';
map_64qam(6,5)='011101';
map_64qam(6,6)='011111';
map_64qam(6,7)='010111';
map_64qam(6,8)='010101';

map_64qam(7,1)='110001';
map_64qam(7,2)='110011';
map_64qam(7,3)='111011';
map_64qam(7,4)='111001';
map_64qam(7,5)='011001';
map_64qam(7,6)='011011';
map_64qam(7,7)='010011';
map_64qam(7,8)='010001';

map_64qam(8,1)='110000';
map_64qam(8,2)='110010';
map_64qam(8,3)='111010';
map_64qam(8,4)='111000';
map_64qam(8,5)='011000';
map_64qam(8,6)='011010';
map_64qam(8,7)='010010';
map_64qam(8,8)='010000';

%%
    %%
       
        
        % Add Noise to Passive Radar
        if (0)
        noisy_sur_channel= awgn(sur_channel,-2,'measured');
        else
            noisy_sur_channel=sur_channel;
        end

        figure,plot(real(data(1:100000)))

        figure,plot(abs(fft(ref_channel(1:100000))))
        
        fs=3.84e9/8/32;
        dt=1/fs;
        delta=[224, 112,56,28]*1e-6;
        dd=delta(4);
        Tu=896e-6;
        D=dd/dt;
        D=413;
        S=Tu/dt;
        Tots=S+D;
        B=S-D;
        P1=zeros(1,100000);

        temp=ref_channel(1:200000);
        %%
        
        for i=1:1:length(temp)-S-D-1
            P1(i)=sum(conj(temp(i:i+D)).*temp(i+S:i+S+D));
        end

        figure
        plot(abs(P1))

        figure
        plot(angle(P1))

        %% getting the indexes of the start of the symbols
        state=0;
        mthreshold=8e5;%0.75*max(abs(P1));
        down_threshold=3e5;
        symbols_starts=[];
        for i=1:1:length(P1)
            if abs(P1(i))>mthreshold
                if state==0
                    current_peak=abs(P1(i));
                    current_index=i;
                    
                elseif abs(P1(i))>current_peak
                    current_peak=abs(P1(i));
                    current_index=i;
                end
                state=1;
                
            else
                if abs(P1(i))<down_threshold && state==1
                    symbols_starts=[symbols_starts,current_index];
                     state=0;
                    current_peak=0;
                end
                % if state==1
                %     symbols_starts=[symbols_starts,current_index];
                % end
                % state=0;
                % current_peak=0;
            end
        end

        %% do frequency aligment 

        freq_avg=mean(angle(P1(symbols_starts)));

        temp2=temp.*exp(+1j*freq_avg);

        P2=zeros(1,100000);
        for i=1:1:length(temp2)-S-D-1
            P2(i)=sum(conj(temp2(i:i+D)).*temp2(i+S:i+S+D));
        end

        freq_avg=mean(angle(P2(symbols_starts)));
        
        %% get the slices of interest 
        i=1;
        XF=fft(temp(symbols_starts(i)+D:symbols_starts(i)+D+S),length(symbols_starts(i)+D:symbols_starts(i)+D+S)-1); 

        figure
        plot(abs(XF))

        %% try this again

temp=ref_channel(1:200000);

for i=1:1:length(temp)-S-D-1
    P1(i)=sum(conj(temp(i:i+D)).*temp(i+S:i+S+D));
end

figure
plot(abs(P1))

% for i=1:1:10
msignal=temp(11982+D+2-1:11982+D+S+1);

XF=fft(msignal);

figure 
scatter(real(XF(1:300)),imag(XF(1:300)))
% end

X_pilots=XF(10:12:end);

figure
plot(abs(X_pilots))

% figure
% for i=1:1:10
XX=linspace(-1-5+7,6819-5+4,length(msignal));
XV=linspace(1,6817,6817);
new_signal=interp1(XX,msignal,XV,"spline");
XF_int=fft(new_signal,6817);
XF_int_temp=fft(new_signal,6817);

figure
hold on
plot(XV,real(new_signal))
plot(XX,real(msignal))

figure 
hold on
scatter(real(XF_int(1:200)),imag(XF_int(1:200)))
scatter(real(XF_int(10)),imag(XF_int(10)),"filled")
scatter(real(XF_int(22)),imag(XF_int(22)),"filled")
scatter(real(XF_int(34)),imag(XF_int(34)),"filled")
scatter(real(XF_int(46)),imag(XF_int(46)),"filled")
scatter(real(XF_int(58)),imag(XF_int(58)),"filled")
scatter(real(XF_int(70)),imag(XF_int(70)),"filled")
scatter(real(XF_int(82)),imag(XF_int(82)),"filled")
scatter(real(XF_int(94)),imag(XF_int(94)),"filled")
scatter(real(XF_int(106)),imag(XF_int(106)),"filled")
scatter(real(XF_int(118)),imag(XF_int(118)),"filled")
scatter(real(XF_int(130)),imag(XF_int(130)),"filled")
scatter(real(XF_int(142)),imag(XF_int(142)),"filled")
scatter(real(XF_int(154)),imag(XF_int(154)),"filled")
scatter(real(XF_int(166)),imag(XF_int(166)),"filled")
scatter(real(XF_int(178)),imag(XF_int(178)),"filled")
scatter(real(XF_int(190)),imag(XF_int(190)),"filled")

% end
XF_int=fftshift(XF_int);
X_pilots_int=XF_int(10:12:end);

figure
hold on
scatter(real(XF_int(1:200)),imag(XF_int(1:200)))
scatter(real(X_pilots_int(1:20)),imag(X_pilots_int(1:20)),"filled")
scatter(real(X_pilots_int(21:40)),imag(X_pilots_int(21:40)),"filled")
scatter(real(X_pilots_int(41:80)),imag(X_pilots_int(41:80)),"filled")
scatter(real(X_pilots_int(81:120)),imag(X_pilots_int(81:120)),"filled")
scatter(real(X_pilots_int(121:150)),imag(X_pilots_int(121:150)),"filled")
scatter(real(X_pilots_int(150:180)),imag(X_pilots_int(150:180)),"filled")

check_prbs=zeros(1,length(X_pilots_int));
check_prbs(imag(X_pilots_int)<0)=1;
figure
plot(abs(X_pilots_int))

prbs_seq = generate_prbs(length(XF_int));

pilot_locations=(4/3)*2*(0.5-prbs_seq)+1j*0;

CS=X_pilots_int./pilot_locations(10:12:end)';

X_pilots_new=X_pilots_int./CS;

figure
hold on
scatter(real(X_pilots_new(1:20)),imag(X_pilots_new(1:20)),"filled")

X_symbols_new=XF_int;
pilots_ind=10:12:length(XF_int);
for i=1:1:6817
    [M,II]=min(abs(i-pilots_ind));
    X_symbols_new(i)=XF_int(i)./CS(II);

end

figure
hold on
scatter(real(XF_int(1:500)),imag(XF_int(1:500)))

figure
hold on
scatter(real(X_symbols_new(1:1500)),imag(X_symbols_new(1:1500)))
figure
hold on
scatter(real(X_symbols_new(1:6817)),imag(X_symbols_new(1:6817)))
scatter(real(X_symbols_new(10)),imag(X_symbols_new(10)),"filled")
scatter(real(X_symbols_new(22)),imag(X_symbols_new(22)),"filled")
scatter(real(X_symbols_new(34)),imag(X_symbols_new(34)),"filled")
scatter(real(X_symbols_new(46)),imag(X_symbols_new(46)),"filled")
scatter(real(X_symbols_new(58)),imag(X_symbols_new(58)),"filled")
scatter(real(X_symbols_new(70)),imag(X_symbols_new(70)),"filled")

scale_factor=6.3636;
sypo=[-7,-5,-3,-1,1,3,5,7]/scale_factor;
i_symbol=zeros(1,6817);
q_symbol=zeros(1,6817);

for i=1:1:6817
   [~,i_symbol(i)]= min(abs(real(X_symbols_new(i))-sypo));
   [~,q_symbol(i)]= min(abs(imag(X_symbols_new(i))-sypo));
end




figure,
hold on
plot(real(new_signal))
plot(imag(new_signal))

temp3=fftshift(XF_int);
% XF_int_temp

figure
hold on
plot(real(temp3))
plot(real(XF_int_temp))

resignal=ifft(ifftshift(XF_int),length(new_signal));

figure
hold on
plot(real(new_signal))
plot(real(resignal))

temp4=ifftshift(XF_int);
temp5=(length(msignal)/6817)*[temp4(1:3409),zeros(1,length(msignal)-6817),temp4(3410:end)];
res2_signal=ifft(temp5);

figure
hold on
plot(real(msignal))
plot(real(res2_signal))

figure
hold on
plot(prbs_seq)
scatter(10:12:length(XF_int),check_prbs)


%%




        % frequency_fix=temp;
        t=1:1:length(temp);
        frequency_fix=exp(1j*0.15/S*t);

        figure,
        hold on
        plot(real(frequency_fix))

        temp=temp'.*frequency_fix;

        for i=1:1:100000
            P1(i)=sum(conj(temp(i:i+D)).*temp(i+S:i+S+D));
        end
        
        % figure
        % hold on
        for j=-2:1:2
            % for i =7:1:11
                figure
                i=8;
                j
                msignal=temp(11982+D:11982+D+S-2);
                XX=linspace(1,6817,length(msignal));
                XV=linspace(1,6817,6817);
                new_signal=interp1(XX,msignal,XV,"spline");
                XF_new=fft(new_signal,6817);            
                scatter(real(XF_new(1:200)),imag(XF_new(1:200)))
                % drawnow
                % pause(1)
                
            % end
        end
        

                msignal=temp(11982+D:11982+D+S-2);
                XX=linspace(1,6817,length(msignal));
                XV=linspace(1,6817,6817);
                new_signal=interp1(XX,msignal,XV,"spline");
                XF_new=fft(new_signal,6817);     
                figure
                hold on
                scatter(real(XF_new(1:200)),imag(XF_new(1:200)))
                scatter(real(XF_new(10)),imag(XF_new(10)),"filled")
                scatter(real(XF_new(22)),imag(XF_new(22)),"filled")
                scatter(real(XF_new(34)),imag(XF_new(34)),"filled")
                scatter(real(XF_new(46)),imag(XF_new(46)),"filled")
                scatter(real(XF_new(58)),imag(XF_new(58)),"filled")
                scatter(real(XF_new(70)),imag(XF_new(70)),"filled")
                scatter(real(XF_new(82)),imag(XF_new(82)),"filled")
                scatter(real(XF_new(94)),imag(XF_new(94)),"filled")


                X_pilots=XF_new(10:12:end);


                figure
                scatter(real(X_pilots(1:200)),imag(X_pilots(1:200)))

                figure
                plot(unwrap(angle(X_pilots)))


        figure
        for i=1:1:10
            XF=fft(temp(11982+D:11982+D+S),length(11982+D:11982+D+S)-10+i);        
            scatter(real(XF(500:700)),imag(XF(500:700)))
            drawnow
            pause(1)

        end

        XF=fft(temp(11982+D:11982+D+S),length(11982+D:11982+D+S)-1); 

        figure
        plot(abs(P1))

        S=13440;

        figure
        hold on
        plot(real(temp(11982:11982+D)))
        plot(real(temp(11982+S:11982+S+D)))
    
        figure
        % for n=1:1:20
        plot((angle(P1)))
        temp2=temp.*exp(-1j*19*pi/100);
        msignal=temp(11982+D:11982+D+S-1);
        XF=fft(temp2(11982+D:11982+D+S),13440);

        XX=linspace(1,6817,length(msignal));
        XV=linspace(1,6817,6817);
        % new_signal=interp1(XX,msignal,XV,"linear");
        new_signal=interp1(XX,msignal,XV,"spline");
        % new_signal=interp1(XX,msignal,XV,"pchip");

        figure
        plot(abs(fft(new_signal)))
        
        XF_new=fft(new_signal);
        figure
        hold on
        scatter(real(XF_new(1:200)),imag(XF_new(1:200))) 

        figure
        hold on
        scatter(real(XF_new(800:1000)),imag(XF_new(800:1000))) 

        figure
        hold on
        scatter(real(XF_new(6000:6200)),imag(XF_new(6000:6200))) 

        figure
        hold on
        scatter(real(XF_new(1201:1400)),imag(XF_new(1201:1400))) 

        figure, hold on
        plot(XX,real(msignal))
        plot(XV,real(new_signal))

        XF_new=fft(new_signal);

        figure
        plot(abs((XF)))

        figure
        hold on
        for i=1:1:floor(length(XF)/12)
            i
            scatter(real(XF(i*12+10)),imag(XF(i*12+10)),"filled")
            drawnow
            pause(0.1)
        end

        % maybe has to be in a different order or like fftshift
        figure
        hold on
        scatter(real(XF(1000:1200)),imag(XF(1000:1200))) %1050 1101 1107 1110
        scatter(real(XF(1050)),imag(XF(1050)),"filled")
        scatter(real(XF(1101)),imag(XF(1101)),"filled")
        scatter(real(XF(1107)),imag(XF(1107)),"filled")
        scatter(real(XF(1110)),imag(XF(1110)),"filled")

        figure
        hold on
        scatter(real(XF(1:200)),imag(XF(1:200)))
        scatter(real(XF(10)),imag(XF(10)),"filled")
        scatter(real(XF(22)),imag(XF(22)),"filled")
        scatter(real(XF(34)),imag(XF(34)),"filled")
        scatter(real(XF(46)),imag(XF(46)),"filled")
        scatter(real(XF(58)),imag(XF(58)),"filled")
        scatter(real(XF(70)),imag(XF(70)),"filled")
        scatter(real(XF(82)),imag(XF(82)),"filled")
        scatter(real(XF(94)),imag(XF(94)),"filled")

        scatter(real(XF(49)),imag(XF(49)),"kx")
        scatter(real(XF(55)),imag(XF(55)),"kx")
        scatter(real(XF(88)),imag(XF(88)),"kx")

        figure
        scatter(real(XF_new(500:700)),imag(XF_new(500:700)))
        scatter(real(XF_new(300:500)),imag(XF_new(300:500)))
        scatter(real(XF_new(100:200)),imag(XF_new(100:200)))
        scatter(real(XF_new(1000:1200)),imag(XF_new(1000:1200)))

        % try resampling the frequency data
        figure
        scatter(real(XF(1:500)),imag(XF(1:500)))
        hold on
        scatter(real(XF(1)),imag(XF(1)),"filled")
        scatter(real(XF(2)),imag(XF(2)),"filled")
        scatter(real(XF(3)),imag(XF(3)),"filled")
        % drawnow 
        % pause(0.5)
        % 
        % end
        figure
        plot(abs((XF)))

        find(real(XF(1:200))<-9000)

        find(real(XF_new(1:200))<-9000)
        
        figure
        scatter(real(XF(1:500)),imag(XF(1:500)))
        scatter(real(XF(2000:2200)),imag(XF(2000:2200)))

       length(XF(abs(XF)>1000))



        % 
        % figure
        % for h=1:1:2000
        %     mstart=h;
        %     mspacing=8192;
        %     XF=fft(ref_channel(mstart:mstart+mspacing));
        %     scatter(real(XF),imag(XF))
        % 
        % end

        % figure
        % plot(abs(fft(ref_channel(1:10000))))
        % 
        %  figure
        % plot(real((ref_channel(:))))
        % 
        % figure,plot(abs(xcorr(ref_channel(1:20000),ref_channel(1:20000))))
        % 
        % 
        % figure,plot(abs(xcorr(ref_channel(1:20000),ref_channel(1:20000),'unbiased')))
        % mshift=2986;
        % mout=sum(abs((ref_channel(1:20000).*conj(ref_channel(1+mshift:20000+mshift)))));
        % 
        % [out] = myxcorr(ref_channel(1:20000),ref_channel(1:20000));
        % 
        % 
        % figure,plot(abs(mout))
        % 
        % figure
        % hold on
        % plot(real(ref_channel(1:20000)))
        % plot(real((ref_channel(1+2987:20000+2987))))
        % 
        % 
        % temp3=zeros(1,40000);
        % for i=1:1:40000
        %     temp3(i)=(sum(ref_channel(i:i+400).*conj(ref_channel(i+2987:i+400+2987))));    
        % end
        % figure,plot(abs(temp3))
        % 
        % figure,plot(angle(temp3))
        % 
        % XF=fft(ref_channel(17994:17994+2987));
        % 
        % figure
        % scatter(real(XF(1:10)),imag(XF(1:10)))
        


        % [ref_matrix_fpga ,self_ambg_matrix_fpga, cc_matrix_fpga] = passive_batch_process_fpga(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
       

        % self_ambg_matrix=self_ambg_matrix_fpga;
        % cc_matrix=cc_matrix_fpga;
        % ref_matrix= ref_matrix_fpga;
        passive.Fs=15e6;
        C=299792458;
        passive.Fc=498e6;
        capture_duration=(length(ref_channel)/passive.Fs);%15;
        ref_channel=ref_channel';
        noisy_sur_channel=noisy_sur_channel';

        [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.min_range,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
        % save(exp_dir + 'passive_matrix','cc_matrix')


        
%% Proccess Passive data into Range-Doppler Slices

        passive.PRF = passive.seg_s; % seg_s
        passive.cpi = 0.5; % coherent proccessing interval (s)
        passive.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
        passive.doppler_window = 'Blackman-Harris';
        passive.dopp_zero_padding = 1;
        passive.dynamic_range = +50;
        passive.range_bin_size = (1/passive.Fs)/passive.range_zero_padding * C/2;
        passive.max_range_m = passive.max_range*passive.range_bin_size;
        passive.min_range_m = passive.min_range*passive.range_bin_size;
        % Create range-Doppler surfaces
        [passive.number_cpi,...
            passive.pulses_per_cpi,...
            passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
            passive.cpi_overlap,...
            passive.dopp_zero_padding,...
            passive.doppler_window);
        % Create Ambiguity surfaces
        [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
            passive.cpi_overlap,...
            passive.dopp_zero_padding,...
            passive.doppler_window);
        
        % Create Active Data Range and Doppler Axis
        
        passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
        passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
        passive.no_range_bins = size(cc_matrix,1);
        passive.range_bins = 1:passive.no_range_bins;
        passive.doppler_bins = passive.pulses_per_cpi*passive.dopp_zero_padding+1;
        passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
        passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
        passive.range_axis = linspace(passive.min_range_m,passive.max_range_m,passive.no_range_bins);
        
        
        % create video of Passive range-Doppler slices        
        video_name = exp_dir + "passive_range-Doppler_" +'simple_running' + ".avi";
        video_title = "Passive Pre-DSI";
        dynamic_range = 80;
        max_range = 100;
        max_doppler = 30;
        frame_rate = 1/(capture_duration/passive.number_cpi);
        createVideo(passive.range_doppler_slices,frame_rate,...
            passive.range_axis,max_range,...
            passive.doppler_axis,max_doppler,...
            dynamic_range,video_name,video_title);
        %
%% Direct Signal Interference Cancellation
% set DSI cancellation parameters

        p = 0.999;        % subtraction parameter - P must be a positive
        % integer less than one to avoid unwanted discontinuities
        % arising from zero values in the range–Doppler surface.
        threshold = 0.005; % cutoff threshold parameter
        max_iterations = 1000; % maximum number of itterations DSI is CLEANed from CAF slice
        number_rbins = size(passive.range_doppler_slices{1},2);
        
        % perform CLEAN based DSI Cancellation
        passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
            passive.self_ambg_slices,...
            ref_matrix,...
            passive.number_cpi,...
            passive.pulses_per_cpi,...
            passive.cpi_stride,...
            passive.dopp_zero_padding,...
            number_rbins,...
            max_iterations,...
            threshold,p,...
            passive.range_axis,passive.doppler_axis);
        
        %   % create video of CLEANed range-Doppler slices
        video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" +'simple_running' + ".avi";
        %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
        video_title = "CLEANed Passive Radar Capture";
        dynamic_range = 80;
        max_range = 100;
        max_doppler = 30;
        frame_rate = 1/(capture_duration/passive.number_cpi);
        createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
            passive.range_axis,max_range,...
            passive.doppler_velocity_axis,max_doppler,...
            dynamic_range,video_name,video_title);        
    end 

toc


% temp_full=passive.CLEANed_range_doppler_slices{45};

% save('image_frame_save','temp_full')

temp=passive.CLEANed_range_doppler_slices{45};

load('image_frame_save')

figure
frame = 20*log10(abs(temp./max(temp(:))));
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,frame, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;

figure
frame_og = 20*log10(abs(temp_full./max(temp_full(:))));
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,frame_og, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;

diff_frame = 20*log10(abs((temp_full-temp)./max(temp_full(:))));
figure
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,diff_frame, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;



