t_signal = 5; % s

aPR = init_aPR(512);
s2 = randn(t_signal * aPR.SampleRate, 2);

[a, b] = butter(5, [100,8000]/(aPR.SampleRate/2), 'bandpass');
s2 = filter(a, b, s2);

output = zeros(t_signal * aPR.SampleRate, 2);
for i = 1:t_signal * aPR.SampleRate / aPR.BufferSize
    input_buff = s2(i * (1:aPR.BufferSize), :);
    [recbuf1] = aPR(input_buff);
    output(i * (1:aPR.BufferSize), :) = recbuf1;
end

t = 0:1/aPR.SampleRate:t_signal - 1 / aPR.SampleRate;
plot(t, output(:, 1))
hold on
plot(t, output(:, 2))

