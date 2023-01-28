function[ysd,OAcoef_out] = GSLC_process_2ChansDelayed(x,OAcoef_in,StAng,Tao,micDist,fs,c)

% [ysd,OAcoef_out] = GSLC_process_2ChansDelayed(x,OAcoef_in,StAng,Tao,micDist,fs,c)
% Computes output of a 2-channel generalized side-lobe cancelling beamformer,
% after receiving channels that have been delay to achieve the desired
% steering angle (Task 2.5)
%
% Inputs:
%
% x:            (Nx2 vector) 2 channel input signal that has had delays applied
%               where N is the length of the input frame
% OAcoef_in:    (Nx1 vector) Overlap-add coefficients from previous call of function. 
%               Should be initialized to a vector of zeros prior to first
%               call.
% StAng:        (scaler) Current steering angle of the beamformer
% Tao:          (scaler) The time difference between delays applied to the
%               two channels (in samples).
% mic_dist:     (scaler) Distance between adjacent array microphones (m)
% fs:           (scaler) Audio sampling frequency (Hz)
% c:            (scaler) speed of sound (m/s)
%
% Outputs:
%
% ysd:         (Nx1 vector) Output of beamformer
% OAcoef_out:  (Nx1 vector) Overlap-add coefficients to be input on next
%               function call.

if size(x,2)>size(x,1)
    x = tranpose(x);
end

xd1 = x(:,1);
xd2 = x(:,2);

if StAng < 0
    ys = 0.5 * (xd1 + xd2); %SDB
    yd = 0.5 * (xd1 - xd2); % cancellation path
else
    ys = 0.5 * (xd1 + xd2); %SDB
    yd = 0.5 * (xd2 - xd1); % cancellation path
end
Lframe = length(x);

YD = zeros(2*Lframe,1); % prealloc fft-space
YD(1:Lframe,:) = yd;    % zeropad signal
YD = fft(YD,2*Lframe);  % calc fft
ta = 2*pi*fs*(1-1/(2*Lframe))*linspace(0,1,2*Lframe); % calc fft frequency vector
WSD = (-1i*sinc((micDist/(c))*ta).*sin(Tao*ta))./(1-sinc((micDist/(c))*ta).*cos(Tao*ta)); % filter calc
WSD(1) = -1i; % account for NaN at 0Hz
ydf = ifft(YD .* WSD',2*Lframe,'symmetric'); % transform to time domain assuming conjugate symmetric signal
%ydf = real(ifft(YD .* WSD',2*frameLen));
olaydf = ydf(1:Lframe) + OAcoef_in;      % overlap-add with last frame
ysd = ys - olaydf;                             %
OAcoef_out = ydf(Lframe+1:end);           % save OLA buffer for next iteration