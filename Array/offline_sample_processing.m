function [playdata] = offline_sample_processing(sample,algo,initdata,initparam,pagesize,fs)
% function [playdata]=offline_sample_processing(sample,algo,initdata,initparam,pagesize,fs)
%
% Offline implementation of the realtime framework used in the class Signal
% Processing for Audio Technology at TUM. In this implementation, a sample
% is used as input and the output is returned as a sample. Other
% implementations:
% realtime_processing_v2-0.m        : real time audio input and output
% realtime_sample_processing_v2-0.m : sample is used as input, output real time
%
% Inputs:
% sample       audio sample to be processed
% algo         plugin object for algorithm implementation (e.g. my_algo)
% initdata     data to be passed to algorithm during 'init' phase
% initparam    cell array with initial values for parameters (optional)
% pagesize     size of input and output buffers used with algorithm
% fs           sample frequency
%
% Outputs:
% playdata     played back data (i.e. output from algorithm)
%
% This function was written for the course Signal Processing for Audio
% Technology. Author: Fritz Menzer, Audio Information Processing, TUM
% v1.1, 15.4.2015: Replaced cell arrays with structs and normal arrays
%                  where suitable; algorithm returns now initial state in
%                  "newstate" rather than "output" (F.M.)
% v2.0, 18.04.2018: Use audio plugins instead of a functions in the same
%                   way as with audioPlayerRecorder built in function of 
%                   Matlab
%                   (N. Kolotzek)

% default values for fs and pagesize
if nargin<6, fs=44100; end
if nargin<5, pagesize=512; end

% determine buffer size and sample rate and compute number of buffers
numpages = ceil(size(sample,1)/pagesize); % process at least the entire sample
if size(sample,1) < pagesize*numpages     % zero-pad sample if necessary
    sample(pagesize*numpages,1) = 0;
end
disp(['starting offline processing (pagesize=' num2str(pagesize) '; fs=' num2str(fs) '; numpages=' num2str(numpages) ')']);

% determine algorithm capabilities
% number of input channels
numinchan = algo.getnuminchan;
if numinchan == -1 % if the algorithm's number of input channels is not fixed
    numinchan = size(sample,2);   % use sample to determine number of channels
elseif numinchan > size(sample,2) % raise error if sample has too few channels
    error(['argument sample must have at least ' num2str(numinchan) ' channels for this algorithm']);
end

% number of output channels
numoutchan = algo.getnumoutchan(numinchan);
if numoutchan == -1 % if the algorithm's number of output channels is not fixed
    numoutchan = 1; % define output to be mono
end
% parameters
pnames = algo.getparamnames;   % get parameter names from algorithm implementation
pranges = algo.getparamranges; % get parameter ranges
param = zeros(size(pranges));      % initialize param array
if nargin >= 4 && numel(initparam) > 0  % initparam is given, use it
    param = initparam;
else                  % otherwise initialize to 0 or closest possible value
    for i = 1:length(pranges)
        param(i) = max(min(0,pranges{i}(2)),pranges{i}(1));
    end
end

% initialize initdata if necessary
if nargin < 3, initdata = []; end

% initialize state
info = struct('initdata',initdata,'fs',fs,'pagesize',pagesize,'numinchan',numinchan,'numoutchan',numoutchan);
algo.initialize(info);

% make a zeros vector for the recorded and played back (=processed) data
playdata = zeros(numpages*pagesize,numoutchan);

% perform offline processing loop
for i = 1:numpages
    % fetch page that was just recorded
    inbuf = sample((i-1)*pagesize+(1:pagesize),1:numinchan);
    % process this buffer
    output = algo.process(inbuf,param);
    % write output to playdata
    playdata((i-1)*pagesize+(1:pagesize),:) = output;
end