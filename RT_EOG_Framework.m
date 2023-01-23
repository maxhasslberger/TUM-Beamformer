% clear the workspace
close all;
clear all
clc;

%% General Intialization Parameters (ADJUST VALUES HERE ONLY) %%%%%%%%%%%%%%
% fixed variables
fs_acq = 500;               % sampling frequency [Hz]
Nacq_buff = 6;              % EOG acquisition buffer size

% variables to be set
runtime = 70;               % set here the runtime of the processing loop in seconds
plot_flag = 1;                % plot on = 1; off = 0;
plot_freq = 5;             % Frequency of plotting [# of acquisition buffers]
is_online = 0;              % online = 1; offline = 0;

calibration_flag = 1;       % allow calibration with specified max_angle
max_angle = 45;             % deg; for calibration (max angle occured)
rec_angle_max = 5;          % deg; max angle to floor signal to 0

% set the paths for including an EOG (from Task 3) and a Calibration file
% (from Task 6)
offline_EOG_file = 'seq_3_data.mat';
calib_file = 'EOG_calib.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create gaze angle estimation file
if ~exist('EOG_calib.mat', 'file')
    calib.trials = 0;
    calib.gradient = 0.0;
    save('EOG_calib.mat', 'calib');
end
    

%% Initialize RT_EOG

% write the variables needed for initialization to state

state.rec_angle_max = rec_angle_max;
state.fs = fs_acq;                
state.Nacq_buff = Nacq_buff;              
state.is_online = is_online;
state.EOG_file = offline_EOG_file;      % give measurement data to RT_EOG.m
state.calib_file = calib_file;

[~,state] = RT_EOG('init',state);

% initialize the buffer outputs with zeros to speed up computation
EOG_Vraw = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);

EOG_V_LP = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 3 LP filtering of raw data
EOG_V_BP = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 3 HP filtering of LP data

EOG_dV = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 4 computing derivative of V_LP
EOG_dVsmooth = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 4 smoothing of derivative of V_LP
EOG_dVthresholded = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 4 applying thresholding to smoothed  derivative of V_LP
EOG_edge_idx = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);     % Task 4 final saccade starts and ends

EOG_delta = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);
EOG_next_delta = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);
EOG_angle = zeros(ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff,1);
% Plotting Routine 1/2: setting up the figure===============================
if plot_flag
    figure
    xlabel('Time (s)')
    %set(gca,'Xlim',[0,runtime],'Ylim',[-1.2 1.2])
    set(gca,'Xlim',[0,runtime])
    hold on
end
%===========================================================================
%% Real-Time Processing Loop
% here you define how many buffers you fetch 
a = tic;
output.delta = 0.0; % debug
init_ignore = 1.5; % s
for fr_idx = 1:ceil(state.fs*runtime/state.Nacq_buff)
%b = tic;
    state.fr_idx = fr_idx;
    %delta_prev = output.V_est; % debug
    [output,state] = RT_EOG('process',state);
    
    EOG_Vraw((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.V_raw;

    EOG_V_LP((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.V_LP;
    EOG_V_BP((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.V_BP;
    
    % Saccade edge detection
    EOG_dV((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.dV;
    EOG_dVsmooth((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.dVsmooth;
    EOG_dVthresholded((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.dVthresholded;
    EOG_edge_idx((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.edge_idx;
    
    % Estimation
    EOG_delta((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.V_est;
    EOG_next_delta((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.next_delta;
    EOG_angle((fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff)) = output.angle;

    %% Plotting Routine 2/2======================================================
    
    if plot_flag
        
        if ~mod(fr_idx,plot_freq)
            idx_offset = (fr_idx-plot_freq)*state.Nacq_buff;
            plot_range = idx_offset +(1:plot_freq*state.Nacq_buff);
            plot(plot_range/state.fs,EOG_Vraw(plot_range),'b')
            

        end        
    end
%  pause_time = (1/fs_acq)-toc(b); %Pause time calculated by buffer fill time - processing and plotting time
%pause_time = 0.001;
 %pause(pause_time);
toc(a)
% if output.delta ~= delta_prev % debug
%     figure
%     plot(EOG_V_LP)
%     hold on;
%     plot(EOG_edge_idx-2)
%     plot(EOG_delta)
%     plot(EOG_next_delta)
%     legend('lp', 'edge', 'delta');
% end
end

%% Stop BIOPAC
biopacAPI(state.is_online,'stopAcquisition')
biopacAPI(state.is_online,'disconnectMPDev')

%% Update estimation file - dynamic update without the need of compute_Vcalib.m
if calibration_flag
    calib = load('EOG_calib.mat').calib;
    calib.gradient = (calib.trials * calib.gradient + max_angle / max(abs(EOG_delta))) / (calib.trials + 1);
    calib.trials = calib.trials + 1;
    save('EOG_calib.mat', 'calib');
end

%% Store signals
data.EOG_V_LP = EOG_V_LP;
data.EOG_V_BP = EOG_V_BP;
data.EOG_edge_idx = EOG_edge_idx-2;
data.EOG_delta = EOG_delta;
data.EOG_angle = EOG_angle;
save('EOG_data.mat', 'data');

%% Plot outputs

t = 0 : 1/fs_acq : length(EOG_Vraw) * 1/fs_acq - 1/fs_acq;

figure
plot(t , EOG_Vraw)
title('Raw EOG signal')
xlabel('Time (s)')

figure
plot(t , EOG_V_LP)
title('LP filtered EOG signal')
xlabel('Time (s)')

figure;
plot(t , EOG_V_BP)
xlabel('Time (s)')
title('BP filtered EOG signal')

figure
plot(t , EOG_dV)
title('Derivative of LP filtered EOG signal')
xlabel('Time (s)')

figure
plot(t , EOG_dVsmooth)
title('Smoothed derivative of LP filtered signal')
xlabel('Time (s)')

%figure
%plot(t , EOG_dVthresholded)
%title('Applied threshold on smoothed derivative of LP filtered signal')
%xlabel('Time (s)')

figure
plot(t , EOG_edge_idx)
title('Saccade starts (+1) and ends (-1)')
xlabel('Time (s)')

%Compound Plots
%Derivative and Smoothed Derivative
figure
plot(t , EOG_dV)
hold on;
plot(t , EOG_dVsmooth)
xlabel('Time (s)')
legend('Derivative of LP filtered signal', 'Smoothed derivative of LP filtered signal');


%BP signal and Saccade Edges
figure;
plot(t , EOG_V_LP)
hold on;
plot(t , EOG_edge_idx-2)
xlabel('Time (s)')
legend ('LP filtered signal', 'Saccade edges')

figure;
plot(t , EOG_V_LP)
hold on;
plot(t , EOG_edge_idx-2)
plot(t , EOG_delta)
plot(t , EOG_next_delta)
% grid on;
xlabel('Time (s)')
legend('V_{LP}' , 'EOG edge idx -2' , 'EOG delta' , 'EOG next delta')

figure;
plot(t , EOG_delta)
xlabel('Time (s)')
title('EOG signal estimation')

figure;
plot(t , EOG_angle)
xlabel('Time (s)')
ylabel('angle (deg)')
title('EOG angle estimation')
