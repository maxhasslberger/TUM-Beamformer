initdata.d = 0.2; % m ; microphone distance
initdata.c = 344; % m/s

algo = DSB;
% algo = DMA;
% algo = DMA_GSLC;

sample = audioread('AIP_Song.wav');
% sample = load('task11.mat').output;

playdata = offline_sample_processing(sample, algo, initdata);

sound(playdata(1:500000, :), 44100)