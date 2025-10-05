classdef QEKF
    %QEKF Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dynamics
        measurements
    end

    methods
        function obj = QEKF(dynamics, measurements)
            %QEKF Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = time_update(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

        function outputArg = measurement_update(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end