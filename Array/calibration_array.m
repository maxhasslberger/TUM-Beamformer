
aPR = init_aPR(512);

output = zeros(2*aPR.SampleRate, 2);
for i = 1:aPR.SampleRate / aPR.BufferSize
    [recbuf1] = aPR(zeros(aPR.BufferSize, 2));
    output(i * (1:aPR.BufferSize), :) = recbuf1;
end

t = 0:1/aPR.SampleRate:2 - 1 / aPR.SampleRate;
% output = output(floor(0.2*aPR.SampleRate):floor(0.6*aPR.SampleRate));
plot(t, output(:, 1))

signal = abs(output(:, 1));
[pks,locs] = findpeaks(signal);

pks_output = zeros(size(output, 1), 1);
pks_output(locs) = pks;

hold on
plot(t, pks_output)

% ampl = mean(pks);
ampl = max(pks);
rms_ampl = ampl / sqrt(2);





