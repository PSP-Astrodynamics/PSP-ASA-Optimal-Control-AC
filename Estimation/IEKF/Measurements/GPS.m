classdef GPS
    %GPS Summary of this class goes here
    %   Left invariant measurement (in global frame)
    
    properties
        sample_rate
        accuracy_horizontal_position
        accuracy_vertical_position
        reference_location
        decay_factor
        gps_sensor
        scale
    end
    
    methods
        function obj = GPS(sample_rate, accuracy_horizontal_position, accuracy_vertical_position, options)
            %GPS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                sample_rate % [Hz]
                accuracy_horizontal_position % [m] noise standard deviation for horizontal position
                accuracy_vertical_position % [m] noise standard deviation for vertical position
                options.reference_location = [40.42721456408179, -86.91406212186251, 187]; % reference latitude, longitude, and altitude [LLA], default is Purdue Bell Tower :)
                options.decay_factor = 0.999 % Matlab default, 0 is white noise, 1 is random walk
                options.seed = 42 % random seed
                options.scale = 1;
            end

            obj.sample_rate = sample_rate;
            obj.accuracy_horizontal_position = accuracy_horizontal_position; % [m]
            obj.accuracy_vertical_position = accuracy_vertical_position; % [m]
            obj.reference_location = options.reference_location;
            obj.decay_factor = options.decay_factor;
            obj.scale = options.scale;

            % Not using velocity, groundspeed, or course outputs from
            % sensor object because they aren't always reliable
            obj.gps_sensor = gpsSensor(ReferenceFrame = "ENU", SampleRate = sample_rate, ReferenceLocation = options.reference_location, HorizontalPositionAccuracy = accuracy_horizontal_position, VerticalPositionAccuracy = accuracy_vertical_position, DecayFactor = options.decay_factor, Seed = options.seed);
        end
        
        function [v] = innovation(obj, X, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Pi = [eye(3), zeros([3, 2])]; % Just for homogenous matrices...
            v = Pi * X.inv.element * [z; 0; 1];
        end

        function [z] = measure(obj, X_true, twist)
            %MEASURE Use GPS sensor model to get noisy measurements of true position
            [position_lla,~,~,~] = obj.gps_sensor(X_true.r' ./ obj.scale', X_true.v' ./ obj.scale');

            position_enu = lla2enu(position_lla, obj.reference_location,'flat')';

            z = position_enu .* obj.scale;
        end

        function [S_inv] = inv_covariance_innovation(obj, X_hat_minus, Sigma_hat_minus, H)
            Pi = [eye(3), zeros([3, 2])];
            S = H * Sigma_hat_minus * H' + obj.N; %Pi * X_hat_minus.inv.element * obj.N * X_hat_minus.inv.element' * Pi';
            S_inv = inv(S);
        end

        function [val] = N(obj)
            % Create measurement noise covariance
            val = diag(([obj.accuracy_horizontal_position, ...
                        obj.accuracy_horizontal_position, ...
                        obj.accuracy_vertical_position] .* obj.scale) .^ 2);
        end

        function [val] = H_l(obj, X)
            val = -[zeros(3, 6), eye(3)];
        end

        function [val] = H_r(obj, X)
            val = H_l * X.inv.Ad;
        end
    end
end

