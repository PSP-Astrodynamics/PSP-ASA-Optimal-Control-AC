classdef Vehicle
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m (1, 1) double % [kg]
        L (1, 1) double % [m]
        I double % [kg m2]
        Len (1,1) double % [m]
        max_gimbal (1, 1) double % [rad]
        min_thrust (1, 1) double % [N]
        max_thrust (1, 1) double % [N]
        Name (1,1) string
    end
    
    methods
        function obj = Vehicle(m, L, Len, max_gimbal, min_thrust, max_thrust, options)
            %VEHICLE Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                m
                L
                Len
                max_gimbal
                min_thrust
                max_thrust
                options.I = Vehicle.estimateMOI(m, Len)
                options.Name =  "vehicle"
            end

            obj.m = m;
            obj.L = L;
            obj.I = options.I;
            obj.Len = Len;
            obj.max_gimbal = deg2rad(max_gimbal);
            obj.min_thrust = min_thrust;
            obj.max_thrust = max_thrust;
            obj.Name = options.Name;
        end

        function save(obj)
            vehicle = obj;
            path = string("Vehicles\" + obj.Name + ".mat");
            save(path, "vehicle");
        end
    end

    methods(Static)
        function I = estimateMOI(m, Len)
            I = (1/12) * m * (Len)^2; 
        end
    end
end

