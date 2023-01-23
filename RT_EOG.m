function[output,state] = RT_EOG(cmd,state)

% RT_EOG
% To be adapted by students of the Ringpraktikum Neuro-Signale course. When
% called, this function should retreive the next availabe EOG acquisition
% buffer from the Biopac MP36 system, process that buffer and estimate the
% angle of the users gaze in real-time. The function is designed to be
% called repeatedly within a real-time loop of EOG_RT_Framework.m.
%
% Input/Output Variables:
%
% cmd:      Input string used to choose which phase of the function is excuted.
%           'init' executes the intialization phase of the function used to
%           compute initial conditions of variables, filter coefficients,
%           etc that should be computed before the real-time loop commences.
%           'process' executes the main processing phases that should be
%           called within the real-time processing loop
%
% state:    Input/Output structure used to store any variables that need to
%           be stored and passed between subsequent calls of RT_EOG. You
%           are free to add any fields to this structure as you see fit,
%           which can be initialized and set within this function. However,
%           the following fields must be included and must be initiliazed
%           OUTSIDE this function.
%
%           state.fs: Sampling rate of EOG acquisition with BIOPAC.
%
%           state.Nacq_buff: Length of EOG acquisition buffer in samples.
%
%           state.is_online: If 1, the function retreives buffers from the
%           BIOPAC system. If 0, function loads in a pre-recorded EOG
%           signal and processes it as if it was being read from the BIOPAC
%
%           state.fr_idx: Tracks how many buffers have been retrieved and
%           processed since the start of execution. Can be used to avoid
%           computing gaze angles while the EOG signal is fluctuating due
%           to system switch-on artefacts.
%
% output:   Output structure comprising an variable you would like to have
%           accessible once the function returns. You can add any fields to
%           this structure you desired, however the following must be
%           included for correct operation of with EOG_RT_Framework
%
%           output.V_raw: the raw EOG data retrieved from the acquistion
%           buffer
%
%           output.V_LP: lowpass filtered version of EOG buffer
%
%           output.V_BP: bandpass version of EOG buffer
%
%           output.edge_idx: vector that is same length as EOG buffer that
%           has contains +1 and -1 on samples where a saccade start and end
%           have respectively been detected, and 0 on others.
%
%           output.ang_est: (scaler) the estimated angle of the user's gaze
%           as of the end of the EOG buffer


if strcmpi(cmd,'process')
    %% Read in next available EOG buffer
%     tic
    [~, output.V_raw] = biopacAPI(state.is_online,'receiveMPData',state.Nacq_buff); % pull eog data frame
%     toc
    %% Task 3: Apply filtering to Acq Buffer
    
        % 1. Low-pass filtering  
        % 2. HP filtering: Apply low cutoff LP and subtract from output of
        % above LP signal

        % 1. LP output = output.V_LP
        [output.V_LP , state.cheby.zf] = filter(state.cheby.b , state.cheby.a , output.V_raw , state.cheby.zf);
            
        % 2. Bandpass out = output.V_BP
        % apply LP Butterworth and subtract from LP filtered raw signal
        [V_butter , state.butter.zf] = filter(state.butter.b , state.butter.a , output.V_LP , state.butter.zf);
        output.V_BP = output.V_LP - V_butter;

    %% Task 4: The saccade edge detection algorithm
        
        % 1. Compute dV/dn
        % 2. Smooth dV/dn
        % 3. Set a detection threshold for the absolute value of dV/dn
        % 4. Perform the subtraction needed to arrive at the edge index output
        %    with values [-1, 0, 1]
        % 5. Write the values that you need for the next iteration of the
        % real-time loop to state
        % 6. exclude first buffers due to oscillations that are not
        % saccades

        % 1.
        % differentiate using always 1 value from previous buffer
        output.dV = diff([state.diffbuff , output.V_LP]); % -> output row vector of length Nacq_buff

        % 2.
        % use FIR moving mean filter to smooth dV
        [output.dVsmooth , state.movmeanfilt.zf] = filter(state.movmeanfilt.b , state.movmeanfilt.a , output.dV , state.movmeanfilt.zf);

        % 3.
        output.dVthresholded = output.dVsmooth;
        output.dVthresholded(abs(output.dVsmooth) < state.saccade.threshold) = 0;
        output.dVthresholded(abs(output.dVsmooth) > state.saccade.threshold) = 1;

        % 4.
        dVsubtract = zeros(1 , state.Nacq_buff);
        dVsubtract(2 : end) = output.dVthresholded(1 : end - 1);
        dVsubtract(1) = state.subbuff;
        output.edge_idx = output.dVthresholded - dVsubtract;

        % 5.
        % diffbuff always keeps end value of previous iteration of V_LP
        state.diffbuff = output.V_LP(end);
        % subbuff always keeps end value of previous iteration of
        % dVthresholded
        state.subbuff = output.dVthresholded(end);

        % 6.
        if state.exclude_i <= state.pausebuff
            output.edge_idx = 0;
            state.exclude_i = state.exclude_i + 1;
        end
        
    %% Task 5: Estimating the potential changes    
    
        % For each sample in the saccade edge sequence
        % state.delta = 0.0;
        state.avg_buff = [state.avg_buff(state.Nacq_buff+1:end), output.V_LP]; % update avg buffer
        
        % 1. Detect the start of a saccade and count how many samples it
        %    takes

        % correct task: detect end of saccade and count samples until next
        % one starts
        onset_id = find(output.edge_idx == 1, 1); % returns index of saccade start in buffer
        offset_id = find(output.edge_idx == -1, 1); % refurns index of saccade end in buffer
        
            % 1: case saccade end is detected and counting hasn't started already:
            % 2: case no saccade start or end in buffer & counting started
            % already:
            % 3: case saccade start:
        if ~isempty(offset_id) && ~state.complete % upcounting
            state.cnt = state.Nacq_buff - offset_id; % start counting with elements after saccade end in buffer
            state.complete = 1; % start further counting
        elseif isempty(onset_id) && isempty(offset_id) && state.complete
            state.cnt = state.cnt + state.Nacq_buff; % proceed counting, add whole buffer to amount of counted samples
        elseif ~isempty(onset_id)
            state.cnt = state.cnt + onset_id; % final cnt, add samples until start is detected
            state.complete = 0; % stop further counting
        end
        
        % 2. ignore sloppy saccades, reset the timer
        if state.cnt < state.t_sloppy * state.fs && ~isempty(onset_id) % = sloppy
            output.edge_idx = zeros(1, state.Nacq_buff); % delete saccade onset because sloppy
            state.cnt = 0; % reset
            state.del_prev_edge = 1; % flag to delete previous saccade end as well
        
        % 3. distinguish between rapid and normal saccades
        elseif state.cnt < state.t_rapid * state.fs && ~isempty(onset_id) && state.first == 0 % = rapid
        % 4. for rapid saccades, the end potential of the current and the 
        %    start potential of the new saccade are the same, and should be
        %    determined via averaging over the number of available samples
        discard_idx = state.cnt - state.t_sloppy * state.fs; % compute which values of buffer to use for averaging as not full buffer can be used
        tmp = mean(state.avg_buff(discard_idx:end)); % perform averaging
        state.delta = state.delta + tmp - state.next_delta; % add potential difference of rapid saccade to overall potential estimate
        state.next_delta = tmp; % starting potential for next saccade
        state.cnt = 0; % reset
        
        
        elseif (state.cnt >= state.t_rapid * state.fs && ~isempty(onset_id)) || (~isempty(onset_id) && state.first == 1) % normal saccade start, computed second as long as not first saccade
        % 5. for normal saccades, the start and end potentials should be
        %    retrieved via averaging over the Tavg previous samples,
        %    respectively
        state.next_delta = mean(state.avg_buff);
        state.avg_once = 1; % flag that start voltage of saccade is known
        state.cnt = 0; % reset
        state.first = 0; % stop second condition of elseif from being true
        
        elseif state.cnt >= state.t_rapid * state.fs && state.avg_once && state.first == 0 % normal saccade end, computed first
            state.delta = state.delta + mean(state.avg_buff) - state.next_delta; % update abolute voltage estimate
            state.avg_once = 0;
        end
        
%         if state.fr_idx > 100
%             output.V_est = state.delta;
%         else
%             output.V_est = 0.0;
%             state.delta = 0.0;
%         end

        output.V_est = state.delta;

        output.next_delta = state.next_delta;
        output.del_prev_edge = state.del_prev_edge;
    
    %% Task 7: Auto recalibration
    if abs(output.V_est) < state.max_floor
        output.V_est = 0.0;
    end
    
    %% Task 6: Convert the potentials to gaze angles
    output.angle = output.V_est * state.EOG_Vgrad;
    
    %% Update All Save Buffers (All Tasks) HERE %%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(cmd,'init')
    %% Task 3: Initialize the filter coefficients & filter state variables for EOG filtering
        
        % 1. calculate the coefficients for the lowpass filter
        % 2. initialize the state of the lowpass filter
        % 3. calculate the coefficients for the highpass filter
        % 4. initialize the state of the highpass filter

        % 1.
        % LP coefficients stored in state.cheby.b and .a
        cheby_at = 60; % attenuation of Chebyshev filter
        cheby_fc = 35; % fc of Chebyshev filter, need to divide by fs/2 for filter input
        [state.cheby.b , state.cheby.a] = cheby2(4 , cheby_at , cheby_fc / (state.fs / 2));    % compute 4th order type 2 Chebyshev filter
        %freqz(state.cheby.b , state.cheby.a,[],state.fs) % for filter visualization
        
        % 2.
        % LP coefficients stored in state.cheby.b and .a
        state.cheby.zf = zeros(1 , max(length(state.cheby.a),length(state.cheby.b))-1);     % initial filter conditions of rRT processing

        % 3.
        % HP coefficients stored in state.butter.b and .a
        butter_fc = 0.16; % fc of Butterworth filter, need to divide by fs/2 for filter input
        [state.butter.b , state.butter.a] = butter(2 , butter_fc / (state.fs / 2));   % compute 2nd order Butterworth filter
        %freqz(state.butter.b , state.butter.a,[],state.fs) % for filter visualization

        % 4. 
        % HP coefficients stored in state.butter.b and .a
        state.butter.zf = zeros(1 , max(length(state.butter.a),length(state.butter.b))-1);     % initial filter conditions of rRT processing

    %% Task 4: Parameters, Coefficients and Save Variables for Saccade Edge Detection
        
        % 1. set the parameters of the moving average filter (length,
        % coefficients, state)
        % 2. set the threshold for dV
        % 3. initialize the values for variables used from previous buffers
        % 4. initialize amount of time for buffers to be excluded from edge
        % detection in the beginning

        % 1.
        % set number of samples used in 1 filtering process
        state.movmeanfilt.filterlength = 13; 
        % compute a and b of FIR filter for filter.m
        state.movmeanfilt.a = 1;
        state.movmeanfilt.b = ones(1 , state.movmeanfilt.filterlength) / state.movmeanfilt.filterlength;
        state.movmeanfilt.zf = zeros(1 , max(length(state.movmeanfilt.a),length(state.movmeanfilt.b))-1);     % initial filter conditions of rRT processing

        % 2.
        state.saccade.threshold = 0.015;
        % Experiment 2 needs a threshold of max 0.001, for Experiment 1 and 3 00.2, for the Clara file 0.015 worked better

        % 3. 
        % scalar to always store last value of previous iteration's
        % output.V_LP for dV computation
        state.diffbuff = zeros(1,1);  
        % scalar to always store last value of previous iteration's
        % output.dVthresholded for delyed subtraction with itself
        state.subbuff = zeros(1,1);

        % 4.
        edgepause = 0.125;   % exclude first 0.25s
        state.pausebuff = round(edgepause * state.fs / state.Nacq_buff); % amount of buffers to be excluded
        state.exclude_i = 1; % initialize counter how many buffers were excluded already
    
    %%  Task 5: Parameters and Save Variables for Saccade Initial and Final Level Estimation
    
        % 1. the number of samples to wait after a saccade end has been
        %    detected to start the potential estimation
        state.t_avg = 0.1; % s      % time required for averaging
        state.t_sloppy = 0.2; % s       % threshold time for sloppy saccades
        state.t_rapid = state.t_avg + state.t_sloppy;       % threshold time for rapid saccaddes
        state.cnt = 0; % samples in last saccade
        state.complete = 0; % 1 if counting is active already, 0 if not
        state.del_prev_edge = 0; % flag to delete previous edge if sloppy saccade was registered
        
        state.slop_Nbuff = ceil(state.t_sloppy * state.fs / state.Nacq_buff); % threshold amount of buffers for sloppy saccade
        state.rapid_Nbuff = ceil(state.t_rapid * state.fs / state.Nacq_buff); % threshold amount of buffers for rapid saccade
        state.avg_Nbuff = ceil(state.t_avg * state.fs / state.Nacq_buff); % amount of buffers for normal averaging
        
        % 2. the minimum number of samples the eye needs to be stationary,
        %    such that no sloppy saccade is detected
        
        
        % 3. the number of samples to average over to estimate the start
        %    and end potential of a saccade
        state.Navg = state.t_avg * state.fs;
        
        % 4. the minimum number of samples the eye must be stationary
        %    before moving, such that no rapid succession saccade is detected
        % 5. the timer triggered once a saccade end has been detected
        % 6. the start and the end potential of the current saccade
        state.first = 1;        % capture voltage at first saccade start

        % 7. the EOG estimate
        state.next_delta = 0;           % potential at saccade start for delta computations
        state.delta = 0;                % overall potential estimate where deltas are added up

        % 8. a buffer of length Tavg (rounded up to the next multiple of
        %    Nacq_buff)
        state.avg_buff = zeros(1, state.Navg);
        
        state.avg_once = 1;  % flag that start voltage of saccade is known
        state.del_prev_edge = 0;
    
    %% Task 6: Calibration gradient
    
        % 1. load voltage-to-gaze-angle calibration file
        % 2. save calibration gradient to state structure
        calib = load(state.calib_file).calib;
        state.EOG_Vgrad = calib.gradient;
    
    %% Task 7: Parameters for auto-recalibration
    
        % 1. the minimum angle estimate that should be rounded to 0
        % 2. the maximum angle estimate that should be rounded to 0
        % 3. calculate the potential estimate for both
        if state.EOG_Vgrad ~= 0.0
            state.max_floor = state.rec_angle_max / state.EOG_Vgrad;
        else
            state.max_floor = 0.0;
        end
    
    %% Initialization Biopac (Do NOT alter)
    if state.is_online
        dirpath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2 Education'; % path of mpdev header file
        biopacAPI(state.is_online,'initMPDevCom',dirpath); % initialize libraries
    else       
        biopacAPI(state.is_online,'initMPDevCom',state.EOG_file);% initialize libraries
    end
    
    biopacAPI(state.is_online, 'help');                     % print available dll functions
    biopacAPI(state.is_online, 'connectMPDev');             % Connect with MP unit
    biopacAPI(state.is_online, 'setSampleRate', state.fs);     % Set sampling rate to 500 Hz
    biopacAPI(state.is_online, 'setAcqChannels',[1 0 0 0]); % Set acquisition channels
    biopacAPI(state.is_online, 'startMPAcqDaemon');         % Start acquisition daemon
    % --> MP device is now ready to record
    biopacAPI(state.is_online, 'startAcquisition');         % called in the end to reduce delaytime upon first buffer pull
    
    output = [];
else
    error('cmd not recognized')
end