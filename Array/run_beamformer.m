clc;
clear;
close all;

realtime = 1; % flag to select between processing recording in real time (1), or reading in audio "sample" instead (0)
time = 1000; % in seconds, duration of realtime processing

% INPUT CALIBRATED MICROPHONE VALUES - REF 94 dBSPL
% g = [6.78e-6 4.78e-6]; % microphone gains, can be loaded to code, but need to correspond with this format
g = [1 1];

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

% switch between online (realtime) and offline data processing
if realtime == 1
    [recdata,playdata] = realtime_processing(aPR,time,algo,initdata,[], [1 2],inputchannels);
else
    sample = audioread('AIP_song.wav');
    [playdata] = realtime_sample_processing(aPR,sample, algo,initdata,[], [1 2]);
end

% eof