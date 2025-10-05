classdef Barometer
    %BAROMETER Summary of this class goes here
    %   Singleton mearsurement ( o _ o ) - only measures part of position
    %   Left invariant singleton measurement (in global frame)
    
    properties
        sample_rate
        accuracy_pressure
        decay_factor
        baro_sensor
        scale
    end
    
    methods
        function obj = Barometer(sample_rate, accuracy_pressure, options)
            %BAROMETER Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                sample_rate % [Hz]
                accuracy_pressure % [Pa] noise standard deviation for pressure
                options.decay_factor = 0.999 % Matlab default, 0 is white noise, 1 is random walk
                options.seed = 42 % random seed
                options.scale = 1
            end

            obj.sample_rate = sample_rate;
            obj.accuracy_pressure = accuracy_pressure; % [Pa]
            obj.decay_factor = options.decay_factor;
            obj.scale = options.scale;

            obj.baro_sensor = barometerSensor(SampleRate = sample_rate, NoiseDensity = accuracy_pressure, DecayFactor = options.decay_factor, Seed = options.seed);
        end
        
        function [v] = innovation(obj, X, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Pi = [eye(3), zeros([3, 2])];
            v = Pi * X.inv.element * [z; 0; 1];
        end

        function [z] = measure(obj, X_true, twist)
            % Measure altitude derived from barometer measurement
            r = X_true.r;
            [~,~,pressure_true,~] = atmoscoesa(r(3) / obj.scale);
            pressure_reading = obj.baro_sensor(pressure_true);
            height_reading = atmospalt(pressure_reading);
            z = height_reading * obj.scale;
        end

        function [S_inv] = inv_covariance_innovation(obj, X_hat_minus, Sigma_hat_minus, H)
            Sigma_tilde = inv(H * Sigma_hat_minus * H');

            R_hat = X_hat_minus.R;

            M = diag([inv(obj.N(X_hat_minus)); 0; 0]);

            S_inv = Sigma_tilde - Sigma_tilde * inv(R_hat * M * R_hat' + Sigma_tilde) * Sigma_tilde;
        end

        function [val] = N(obj, X)
            % Create measurement noise covariance
            step_size = 1e-2; % [m]

            r = X.r;

            % Numerically take derivative of height w.r.t. pressure
            h_expected = r(3) ./ obj.scale;
            h_stepped = h_expected + step_size;
            [~,~,p_expected,~] = atmoscoesa(h_expected);
            [~,~,p_stepped,~] = atmoscoesa(h_stepped);
            
            partial_h_partial_p = step_size / (p_stepped - p_expected) .* obj.scale;

            % Convert pressure variance to height variance
            val = (partial_h_partial_p * obj.accuracy_pressure) ^ 2;
        end

        function [val] = H_l(obj, X)
            val = -[zeros(3, 6), eye(3)];
        end

        function [val] = H_r(obj, X)
            val = H_l * X.inv.Ad;
        end
    end
end

