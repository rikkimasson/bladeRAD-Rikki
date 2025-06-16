% load signal and split ref and sur
clear 
clc
close all
addpath('../generic_scripts/matlab',...
        '../generic_scripts',...
        '../generic_scripts/ref_signals/') % path to generic functions


addpath('../../RADCOM/Binary_Generation/chirplet/Chirplet_phy/')
% exp_dir= "~/Documents/bladerad_data/Captures/5/"; 
exp_dir= "D:/Captures/11/"; 
Experiment_ID = 0011;       % Experiment Name

graphs=1;

file_location = exp_dir + 'active_' + Experiment_ID;
filename= file_location + ".sc16q11";
% [ref_channel, sur_channel]  = load_passive_data(file_location);
[ signal, signal_i, signal_q ] = load_sc16q11(filename);
% Plot time domain signals
aaa=23445;


% figure
% plot(real(Y_out))
% axis([129400,129800,-inf, inf])
% xlabel('Time (\mus)','FontSize',15)
% ylabel('Frequency (MHz)','FontSize',15, 'Interpreter', 'tex')
% box on
% set(gca,'FontSize',15)
% set(gcf,'color','w');
% set(gca,'linewidth',2)
% % legend('Envelope pulse', 'Peak pulse','1st min pulse')
% % legend('boxoff')
% axis([-inf inf, -inf ,50])
% grid on
% set(gca,'XMinorTick','on','YMinorTick','on')
% hold off
% set(gcf,'units','inches')
% set(gca,'linewidth',2)
% width=8;
% height=5;
% set(gcf,'defaultFigurePaperPositionMode','manual')
% set(gcf, 'Position',[0 0 width height],'PaperSize',[width,height],'PaperPositionMode','auto','InvertHardcopy','off','Renderer','painters');
% exportgraphics(gca,'fm_fmcw_inst_freq.pdf')
if graphs
    figure
    hold on
    plot(signal_i(1:10:length(signal_i)))
end

signal=signal(2516210:3061190);

if graphs
    figure
    hold on
    plot(real(signal))
    plot(imag(signal))
    xlabel('Time samples','FontSize',15)
    ylabel('Normalized Magnitude','FontSize',15, 'Interpreter', 'tex')
    box on
    set(gca,'FontSize',15)
    set(gcf,'color','w');
    set(gca,'linewidth',2)
    % legend('Envelope pulse', 'Peak pulse','1st min pulse')
    % legend('boxoff')
    % axis([-inf inf, -inf ,50])
    grid on
    set(gca,'XMinorTick','on','YMinorTick','on')
end

if graphs
    [output] = display_frequency_material(signal,1063);
end
BW=20e6;
Tc=500e-6;
dt=1/25e6;
SUB_BANDS=5;
fc=0;
FFT_SIZE=512;
SNR=40;
sampled_bw=25e6;

% test_signal=awgn(signal,SNR,'measured');
test_signal=signal;

%% turning stuff to binary file stuff
addpath("C:\Users\rw\Documents\RADCOM\Binary_Generation\Beamformer_testing")

% II=real(test_signal);
% QQ=imag(test_signal);
% 
% 
% IQIQ=zeros(1,2*length(II));
% IQIQ(1:2:length(IQIQ))=II;
% IQIQ(2:2:length(IQIQ))=QQ;
% decimated_signal = decimate(test_signal,4);

% write_binary_to_text('input_comms_signal',test_signal)








%% end of producing binary file stuff

[chirp_order_total,chirp_up_down] = chirplet_phy(test_signal,BW,Tc,dt,fc,SUB_BANDS,FFT_SIZE,graphs,sampled_bw);
% chirp_order_total=chirp_order_total-1;
% total_bits=[chirp_order_total,chirp_up_down];
% temp=binaryVectorToHex(total_bits)
load('temporay_file.mat')
input_sequence=temp;
packet_length=120;
[outputArg1,outputArg2] = return_packets_from_sequence(input_sequence,packet_length);


% packet_length=120;
% pos_in_packet=0;
% look_for_start_sequence=1;
% start_sequence=0;
% enter_packet=0;
% packet={};
% number_packets=0;
% for i=1:1:length(temp)
%     if look_for_start_sequence
%         if temp{i}=='4'
%             start_sequence=start_sequence+1;
%             if start_sequence>7
%                 enter_packet=1;
%                 look_for_start_sequence=0;
%                 start_sequence=0;
%             end
%         else
%             start_sequence=0;
%         end
%     else
%         if enter_packet
%             if pos_in_packet<packet_length
%                 packet=[packet,temp{i}];
%                 pos_in_packet=pos_in_packet+1;
%             else
%                 % end of packet
%                 packet=[packet,temp{i}];
%                 enter_packet=0;
%                 look_for_start_sequence=1;
%                 pos_in_packet=0;
%                 writecell(packet,sprintf("packet%d.txt",number_packets),'Delimiter',',');
%                 number_packets=number_packets+1;
%                 packet={};
%             end
%         end
%     end
% 
% end