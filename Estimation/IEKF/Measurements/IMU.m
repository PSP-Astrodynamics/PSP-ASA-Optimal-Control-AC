classdef IMU
    %IMU Summary of this class goes here
    %   Right invariant measurement (in local frame)
    
    properties
        sample_rate
        accel_parameters
        gyro_parameters
        mag_parameters
        imu_sensor
        magnetic_field
        scale
    end
    
    methods
        function obj = IMU(sample_rate, accel_parameters, gyro_parameters, mag_parameters, options)
            %IMU Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                sample_rate % [Hz]
                accel_parameters % accelparams object describing accelerometer
                gyro_parameters % gyroparams object describing gyroscope
                mag_parameters % magparams object describing magnetometer
                options.seed = 42 % random seed
                options.magnetic_field = [27.555, -2.4169, -16.0849] % [micro Tesla] Matlab default - could use [XYZ,H,D,I,F] = wrldmagm(height,latitude,longitude,decimalYear) where XYZ is in NED
                options.scale = 1
            end

            obj.sample_rate = sample_rate;
            obj.accel_parameters = accel_parameters; % [m]
            obj.gyro_parameters = gyro_parameters; % [m]
            obj.mag_parameters = mag_parameters;
            obj.magnetic_field = options.magnetic_field; % [micro Tesla]
            obj.scale = options.scale;

            obj.imu_sensor = imuSensor('accel-gyro-mag', ReferenceFrame = "ENU", SampleRate = sample_rate, Accelerometer = accel_parameters, Gyroscope = gyro_parameters, Magnetometer = mag_parameters, MagneticField = options.magnetic_field, Seed = options.seed);
        end

        function [z] = measure(obj, X_true, twist)
            %MEASURE Use IMU sensor model to get noisy measurements of true
            % acceleration, angular velocity, and local magnetic field
            [accel_reading, gyro_reading, mag_reading] = obj.imu_sensor(-(twist.a ./ obj.scale(1:3))' - [0, 0, 9.81], (twist.w ./ obj.scale(4:6))', X_true.R'); % Needs orientation of IMU w.r.t. inertial frame

            accel_reading = accel_reading' + X_true.R' * [0; 0; 9.81];

            z = [accel_reading; gyro_reading'; mag_reading'] .* obj.scale;
        end

        % Bad version
        % function [z] = measure(obj, X_true, twist)
        %     %MEASURE Use IMU sensor model to get noisy measurements of true
        %     % acceleration, angular velocity, and local magnetic field
        %     [accel_reading, gyro_reading, mag_reading] = obj.imu_sensor((twist.a ./ obj.scale(1:3))', (twist.w ./ obj.scale(4:6))', X_true.R); % Needs orientation of IMU w.r.t. inertial frame
        % 
        %     z = [accel_reading'; gyro_reading'; mag_reading'] .* obj.scale;
        % end

        % Need to define how IMU would be used for measurement update
        % Actually need to add angular velocity and acceleration into state
        function [v] = innovation(obj, X, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Pi = [eye(3), zeros([3, 2])]; % Just for homogenous matrices...
            v = Pi * X.inv * [z; 0; 1];
        end

        function [S_inv] = inv_covariance_innovation(Sigma_hat_minus, H)
            %S = H * Sigma_hat_minus * H' + N;
            %S_inv = inv(S);
        end

        function [val] = H_l(obj, X)
            %val = -[zeros(3, 6), eye(3)];
        end

        function [val] = H_r(obj, X)
            %val = H_l * X.inv.Ad;
        end


        % Magnetometer Update
        function [v] = mag_innovation(obj, X, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Pi = [eye(3), zeros([3, 2])]; % Just for homogenous matrices...
            v = Pi * X.element * [z; 0; 1] - obj.magnetic_field'; % Multiplying by X rotates into world frame
        end
        
        function [val] = mag_N(obj)
            val = diag(obj.mag_parameters.NoiseDensity .^ 2);
        end

        function [S_inv] = mag_inv_covariance_innovation(obj, X_hat_minus, Sigma_hat_minus, H)
            S = H * Sigma_hat_minus * H' + obj.mag_N; % PROBABLY NEED TO MULTIPLY NOISE DENSITY BY sqrt(obj.sample_rate)!!!!
            S_inv = inv(S);
        end

        function [val] = mag_H_l(obj, X)
            val = obj.mag_H_r * X.inv.Ad();
        end

        function [val] = mag_H_r(obj, X)
            val = [-skew(obj.magnetic_field), zeros(3, 6)]; % make sure correct, check algebra
        end
    end
end

