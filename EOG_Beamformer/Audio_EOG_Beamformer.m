%%%%% EOG

% clear the workspace
close all;
clear all
clc;

%% General Intialization Parameters (ADJUST VALUES HERE ONLY) %%%%%%%%%%%%%%
% fixed variables
fs_acq = 500;               % sampling frequency [Hz]
Nacq_buff = 6;              % EOG acquisition buffer size

% variables to be set
runtime = 1000;               % set here the runtime of the processing loop in seconds
plot_flag = 1;                % plot on = 1; off = 0;
plot_freq = 5;             % Frequency of plotting [# of acquisition buffers]
is_online = 0;              % online = 1; offline = 0;

calibration_flag = 1;       % allow calibration with specified max_angle
max_angle = 45;             % deg; for calibration (max angle occured)
rec_angle_max = 5;          % deg; max angle to floor signal to 0

% set the paths for including an EOG (from Task 3) and a Calibration file
% (from Task 6)
% offline_EOG_file = 'seq_3_data.mat';
% calib_file = 'EOG_calib.mat';
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

%%%%%% Initialize Beamformer

realtime = 1; % flag to select between processing recording in real time (1), or reading in audio "sample" instead (0)
time = 1000; % in seconds, duration of realtime processing
g = [6.78e-6 4.78e-6]; % microphone gains, can be loaded to code, but need to correspond with this format
algo = DMA_GSLC; % to select beamformer implementation: delay_and_sum (Task 1), differential_microphone_array (Task 2), or generalized_sidelobe_canceller (Task 3)

initdata.BufferSize = 512; % initialize buffersize

% set beamformer parameters
initdata.Nmic = 2;
initdata.micDist = 0.2;
initdata.g = g;

initdata.d = initdata.micDist; % m ; microphone distance
initdata.c = 344; % m/s

if length(g) ~= initdata.Nmic
    error('Number of supplied mic gains needs to match number of mics')
elseif size(g(:),1) == size(g,1)
    g = g.';
end

% this is just to remove hardcoded values from the DMA (not needed for the other algos)
inputchannels = 1:initdata.Nmic;
initdata.sum_sign = ones(size(inputchannels));
initdata.sum_sign(1:round(length(initdata.sum_sign)/2)) = -1;

[aPR] = init_aPR(initdata.BufferSize); % initialize audioplayerrecorder (M-Audio Drivers need to be installed)