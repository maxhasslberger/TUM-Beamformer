classdef dB_panning% < audioPlugin
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
%             proc_buf = input.*param(1);

            % Apply gain
%             r_rms = rms(input);
            gain = [10^(-param(1)/20), 10^(param(1)/20)];
            proc_buf = gain .* input; 
            
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
            plugin.state.initdata = initdata.initdata;
            plugin.state.fs = initdata.fs;
            plugin.state.BufferSize = initdata.pagesize;
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
%             output = -1; % arbitrary number of input channels
            output = 2; % arbitrary number of input channels
        end        
        % get parameter names
        function output = getparamnames(plugin)
            % return a cell array with parameter names
            % Provides labels for the slider bars controlling each parameter -
            % One slider bar for each parameter
            % Should have same number of elements as 'getparamranges'
            output = {'gain'};
        end
        % get parameter range
        function output = getparamranges(plugin)
            % return a cell array with parameter ranges
            % Set the min and max values of the slider bar of each parameter
            % Should have same number of elements as 'getparamnames'
            output = {[-5,5]};
        end
    end
end