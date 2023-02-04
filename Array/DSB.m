classdef DSB < handle
    properties
        % initialize state structure
        state = struct();
        % initialize output
        output = [];
    end
    properties (Access = private)
        % insert your private parameters here. These parameters are not
        % visible for the user.
        % No changes needed
    end
    methods
        % processing algorithm
        function proc_buf = process(plugin, input, param)
            % calculate output from input, state, and param
            % to have access to saved additional parameters (e.g. var) in 
            % state type in var = plugin.state.x;
            
            % Preamplify input signal
            input = input ./ plugin.state.id.g;

            % Prefilter input signal
            [input , plugin.state.zf] = filter(plugin.state.b , plugin.state.a , input , plugin.state.zf);

            % Compute delay samples + distant microphone
            try
            plugin.state.angle = load('Angle_transfer.mat').eog_angle;
            end
            plugin.state.angle
            delay_n = round(plugin.state.fs * plugin.state.delay(plugin.state.angle));


            mic_idx = delay_n > 0.0; % Delay mic 2, if delay positive
            delay_n = abs(delay_n);
            
            % Compute output
            input_ext = [plugin.state.save_buf; input]; % update buffer
            sz = size(input, 1);
            proc_buf = (input(:, 2-mic_idx) + input_ext(end-delay_n - sz+1 : end-delay_n, mic_idx+1)) / 2;
            
            % update save buffer and copy signal to second output channel
            plugin.state.save_buf = input_ext(1 + sz : end, :);
            proc_buf = repmat(proc_buf, 1, size(input, 2));
            
            
            
            
            % update state
            % e.g. plugin.state.x = plugin.state.x + 1;
            
            % write output to outbuf to be available outside the function
            plugin.output = proc_buf;
        end
        % initialization routine
        function initialize(plugin,initdata)
            % generate a default state and return it
            % e.g. plugin.state.BufferSize = initdata.pagesize;
            % THINK OF IMPLEMENTING THE DELAYLINE FORMULAR EQ. 1 P. 31
            plugin.state.id = initdata.initdata;
            plugin.state.fs = initdata.fs;
            plugin.state.BufferSize = initdata.pagesize;
            plugin.output = zeros(plugin.state.BufferSize, 2);
            
            % based on 2 mics
            plugin.state.delay = @(phi)plugin.state.id.d / plugin.state.id.c * sin(phi/180*pi);
            plugin.state.save_buf = zeros(ceil(initdata.fs * plugin.state.id.d / plugin.state.id.c), 2);
            
            plugin.state.angle = 0;
            
            % Init filter param
            butter_fc = 1500;
            [plugin.state.b , plugin.state.a] = butter(2 , butter_fc / (initdata.fs / 2), 'low');
%             butter_fc = 4000;
%             [plugin.state.b , plugin.state.a] = butter(2 , butter_fc / (initdata.fs / 2), 'high');
%             butter_fc = [1500, 4000];
%             [plugin.state.b , plugin.state.a] = butter(2 , butter_fc / (initdata.fs / 2), 'bandpass');
            
            plugin.state.zf = zeros(1 , max(length(plugin.state.a),length(plugin.state.b))-1);
        end
        % get number of input channels
        function output = getnuminchan(plugin)
            % return number of input channels
            % No changes needed
            output = -1; % arbitrary number of input channels
        end
        % get number of output channels
        function output = getnumoutchan(plugin,input)
            % return number of input channels
            % No changes needed
            % input here number of channels
            output = -1; % arbitrary number of input channels
%             output = 2; % arbitrary number of input channels
        end        
        % get parameter names
        function output = getparamnames(plugin)
            % return a cell array with parameter names
            % Provides labels for the slider bars controlling each parameter -
            % One slider bar for each parameter
            % Should have same number of elements as 'getparamranges'
            output = {'gain', 'angle'};
        end
        % get parameter range
        function output = getparamranges(plugin)
            % return a cell array with parameter ranges
            % Set the min and max values of the slider bar of each parameter
            % Should have same number of elements as 'getparamnames'
            output = {[-5,5], [-90,90]};
        end
    end
end