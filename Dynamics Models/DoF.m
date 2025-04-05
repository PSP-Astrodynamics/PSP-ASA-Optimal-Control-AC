classdef DoF 
    %DOF Properties of dynamics models
    %   This enumerates all of the current dynamics models and stores their
    %   characteristics for ease of access.
    
    properties
        nx % Number of states
        nu % Number of controls
        ir % Indices of the position vector
        iv % Indices of the velocity vector
        itheta % Indices of the orientation vector
        iw % Indices of the angular velocity vector
        ithrust % Indices of the thrust vector
        constraints % Cell array of the additional constraints for the model
    end
    
    methods
        function obj = DoF(nx, nu, ir, iv, itheta, iw, ithrust)
            %DOF Construct an instance of this class
            %   Detailed explanation goes here
            obj.nx = nx;
            obj.nu = nu;
            obj.ir = ir;
            obj.iv = iv;
            obj.itheta = itheta;
            obj.iw = iw;
            obj.ithrust = ithrust;

            if nx == 6
                obj.constraints = {@(x, u, vehicle) u(:, 2).^2 <= vehicle.max_gimbal^2};
            else
                obj.constraints = {@(x, u, vehicle) u(:, ithrust(1)) > cos(vehicle.max_gimbal) ./ vecnorm(u, 2, 2)};
            end
        end
    end
    enumeration 
        planar3DoF (6, 2, 1:2, 3:4, 5, 6, 1)
        euler5DoF (11, 3, 1:3, 4:6, 7:9, 10:11, 1:3)
        quat5DoF (12, 3, 1:3, 4:6, 7:10, 11:12, 1:3)
        euler6DoF (12, 4, 1:3, 4:6, 7:9, 10:12, 1:3)
        quat6DoF (13, 4, 1:3, 4:6, 7:10, 11:13, 1:3)
    end
end