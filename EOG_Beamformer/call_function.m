initdata.d = 0.2; % m ; microphone distance
initdata.c = 344; % m/s

% algo = DSB;
% algo = DMA;
 algo = DMA_GSLC;
% algo = dB_panning;

% sample = audioread('AIP_Song.wav');
% sample = load('task11.mat').output;
 sample = load('recording.mat').recording_buf;

 
% playdata = offline_sample_processing(sample, algo, initdata);
playdata = realtime_processing(aPR, 10, algo, initdata);

% sound(playdata(1:500000, :), 44100)
sound(playdata(:, :), 44100)