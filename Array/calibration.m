
aPR = init_aPR(512);
% noise = randn(aPR.SampleRate, 2);
% noise = noise / max(abs(noise(:)));
% [~] = aPR(noise);

playbuf = zeros(aPR.BufferSize, 2);
% playbuf(1, :) = 1;
% t = 0:1/aPR.SampleRate:(aPR.BufferSize-1) / aPR.SampleRate;
% playbuf = sin(2*pi*(t-aPR.BufferSize/2/aPR.SampleRate)*5000) ./ (2*pi*(t-aPR.BufferSize/2/aPR.SampleRate) * 5000);
% playbuf(2,:) = sin(2*pi*(t-aPR.BufferSize/2/aPR.SampleRate)*5000) ./ (2*pi*(t-aPR.BufferSize/2/aPR.SampleRate) * 5000);
% playbuf = playbuf';
% plot(playbuf)

input = zeros(2*aPR.SampleRate, 2);
% input(1:aPR.BufferSize:end, :) = 1;
output = zeros(2*aPR.SampleRate, 2);

for i = 1:aPR.SampleRate / aPR.BufferSize
%     playbuf = randn(aPR.BufferSize, 2);
%     playbuf = playbuf / max(abs(playbuf(:)));
    if i ~= 2
        playbuf(1, :) = 0;
    else
        playbuf(1, :) = 1;
    end
    [recbuf] = aPR(playbuf);
    input(i * (1:aPR.BufferSize)) = playbuf(:, 1);
    output(i * (1:aPR.BufferSize)) = recbuf;
end
    
t = 0:1/aPR.SampleRate:2 - 1 / aPR.SampleRate;
plot(t, input(:, 1))
hold on;
plot(t, output(:, 1))
% plot(output(:, 1) / max(abs(output(:, 1))))
release(aPR)
