function [playdata]=realtime_sample_processing(aPR,sample,algo,initdata,initparam,playChanList)
% function [playdata]=realtime_sample_processing(aPR,sample,algo,initdata,initparam,playChanList)
%
% Realtime framework used in the class Signal Processing for Audio
% Technology at TUM. In this implementation, a sample is used as input and
% output is played in realtime using the Audio System Toolbox. Other 
% implementations:
% realtime_processing.m        : real time audio input and output
% offline_processing.m         : input and output are samples
%
% Inputs:
% aPR          audioPlayerRecorder (Audio System Object)
% sample       audio sample to be processed
% algo         plugin object (e.g. a filter plugin)
% initdata     data to be passed to algorithm during 'init' phase
% initparam    cell array with initial values for parameters (optional)
% playChanList vector containing the channel numbers for playback
%
% Outputs:
% playdata     played back data (i.e. output from algorithm)
%
% This function was written for the course Signal Processing for Audio
% Technology. Author: Fritz Menzer, Audio Information Processing, TUM
% v1.1, 15.4.2015: Replaced cell arrays with structs and normal arrays
%                  where suitable; algorithm returns now initial state in
%                  "newstate" rather than "output" (F.M.)
% v1.2, 12.4.2018: Use MATLAB built-in toolbox instead of playrec.
%                  Now the Audio System Toolbox from MATLAB (since R2017a)
%                  is used (audioPlayerRecorder & audioDeviceReader)
%                  (N. Kolotzek)
% v1.3, 18.04.2018 Instead of algo, we're now using audio Plugins
%                  (N. Kolotzek)

pagesize = aPR.BufferSize;
fs = aPR.SampleRate;
numpages=ceil(size(sample,1)/pagesize); % process at least the entire sample
if size(sample,1)<pagesize*numpages     % zero-pad sample if necessary
    sample(pagesize*numpages,1)=0;
end
disp(['starting realtime processing (pagesize=' num2str(pagesize) '; fs=' num2str(fs) '; numpages=' num2str(numpages) ')']);

numinchan=algo.getnuminchan;
if numinchan == -1 % if the algorithm's number of input channels is not fixed
    numinchan = size(sample,2);   % use sample to determine number of channels
elseif numinchan>size(sample,2) % raise error if sample has too few channels
    error(['argument sample must have at least ' num2str(numinchan) ' channels for this algorithm']);
end
% number of output channels
numoutchan=algo.getnumoutchan(numinchan);
if numoutchan==-1 % if the algorithm's number of output channels is not fixed
    if nargin>=6  % if output channels were provided, use them to determine number of channels.
        numoutchan = length(playChanList);
        if numinchan == 1
            out_dist = ones(1,length(playChanList));
        end
        % change the audioPlayerRecorder object to provided playChanList
        release(aPR);
        aPR.PlayerChannelMapping = playChanList;
    else% otherwise use the initialized setting of the object (default from init_aPR: 2 playback channels at channels [1 2])
        numoutchan = length(aPR.PlayerChannelMapping);
        playChanList = aPR.PlayerChannelMapping;
        if numinchan == 1
            out_dist = ones(1,length(playChanList));
        end
    end
else             % otherwise, if the number of output channels is fixed
    if nargin>=6 % raise an error if user provided a different number of channels
        if numoutchan~=length(playChanList)
            if numoutchan == 1
                out_dist = ones(1,length(playChanList));
            else
                error(['argument playChanList must have ' num2str(numoutchan) ' elements for this algorithm']);
            end
        end
    else         % if no channels were provided, use channels 1 to numoutchan.
        playChanList = 1:numoutchan;
        % change the audioPlayerRecorder object to playChanList
        release(aPR);
        aPR.PlayerChannelMapping = playChanList;
    end
end

% parameters
pnames=algo.getparamnames;   % get parameter names from algorithm implementation
pranges=algo.getparamranges; % get parameter ranges
param=zeros(size(pranges));      % initialize param array
if nargin>=5 && numel(initparam)>0,  % initparam is given, use it
    param=initparam;
else                  % otherwise initialize to 0 or closest possible value
    for i=1:length(pranges)
        param(i)=max(min(0,pranges{i}(2)),pranges{i}(1));
    end
end

% initialize initdata if necessary
if nargin<4, initdata=[]; end

% initialize state
info=struct('initdata',initdata,'fs',fs,'pagesize',pagesize,'numinchan',numinchan,'numoutchan',numoutchan);
algo.initialize(info);

% initialize GUI
if ~isempty(pnames)
    % close figure used previously for this algorithm and set up a new one.
    close(findobj('type','figure','name',['realtime_processing: ' class(algo)]))
    fignum=figure('name',['realtime_processing: ' class(algo)]); % give it a name so it can be found again
    % create vectors for slider and label handles
    paramlabels=zeros(size(pnames));
    sliders=zeros(size(pnames));
    valuelabels=zeros(size(pnames));
    % create the sliders and labels (names and values)
    for i=1:length(pnames)
        paramlabels(i) = uicontrol('Style', 'text', 'String',pnames{i},'Position', [-45+i*60,390,50,20]);
        sliders(i) = uicontrol('Style', 'slider', 'Min',pranges{i}(1),'Max',pranges{i}(2),'Value',param(i),...
            'Position', [-30+i*60,40,20,340]);
        valuelabels(i) = uicontrol('Style', 'text', 'String',num2str(param(i)),'Position', [-45+i*60,10,50,20]);
    end
    % set size of figure window
    pos=get(fignum,'Position');
    pos(3)=length(pnames)*60+20; % width
    pos(4)=420;                  % height
    set(fignum,'Position',pos);
    % make sliders and labels adapt to changing window size
    set(paramlabels,'Units','normalized');
    set(sliders,'Units','normalized');
    set(valuelabels,'Units','normalized');
    drawnow limitrate
end

% make a zeros vector for the played back (=processed) data if necessary
if nargout>=1
    playdata=zeros(numpages*pagesize,numoutchan);
end

% perform real time processing loop
for i=1:numpages
    % update parameters
    drawnow limitrate % update parameter values
    for j=1:length(param)
        param(j)=get(sliders(j),'Value'); % read value from slider
        set(valuelabels(j),'String',num2str(param(j))); % update label
    end
    drawnow limitrate % needs to be called again to draw the updated labels

    % fetch page that was just recorded
    inbuf=sample((i-1)*pagesize+(1:pagesize),1:numinchan);
    
    % audio stream loop
    if i==1
        % in the first step, record and play silence in the next buffer
        [~] = aPR(zeros(pagesize,2));
    else
        % in the further steps:
        % use 'process'-method of plugin to process audio data while
        % recording.
        if size(sample,2)==1
            % if sample is mono play back stereo
            [~] = aPR(algo.process(inbuf,param)*ones(size(playChanList)));
            % Write out_buf of plugin to playdata
            playdata(1+(i-1)*pagesize:i*pagesize,1:numoutchan) = algo.output*ones(size(playChanList));
        elseif size(sample,2)==2
            % if sample is stereo play back stereo
            [~] = aPR(algo.process(inbuf,param));
            % Write out_buf of plugin to playdata
            playdata(1+(i-1)*pagesize:i*pagesize,1:numoutchan) = algo.output;
        else
            error('Input signal is not a mono nor a stereo signal!');
            
        end
    end   
end