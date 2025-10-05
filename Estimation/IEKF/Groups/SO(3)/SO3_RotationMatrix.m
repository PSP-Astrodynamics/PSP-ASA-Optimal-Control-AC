classdef SO3_RotationMatrix < SO3
    %SO3_ROTATIONMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "SO(3) Rotation Matrix"
        size = 9
        element
        identity = eye(3)
    end
    
    methods
        function obj = SO3_RotationMatrix(R)
            %SO3_ROTATIONMATRIX Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = R;
        end

        function val = R(obj)
            val = obj.element;
        end
        
        function composition = compose(X, Y)
            % composition = X.element * Y.element - just matrix multiplication
            composition = SO3_RotationMatrix(X.R * Y.R);
        end
        function inverse = inv(X)
            inverse = SO3_RotationMatrix(X.R');
        end

        function val = constraint(X)
            val = X.R' * X.R;
        end
        function val = act(X, p)
            val = X.R * p;
        end

        function tau = vee(G, tau_hat)
            theta = [tau_hat(3, 2); tau_hat(1, 3); tau_hat(2, 1)];
            tau = theta;
        end
        function tau_hat = hat(G, tau)
            theta = tau(1:3);
            tau_hat = skew(theta);
        end

        function tau = Log(X)
            theta = acos((trace(X.R) - 1) / 2);
            tau = theta * X.vee(X.R - X.R') / (2 * sin(theta));
        end

        function X = cayley(G, tau)
            theta = norm(tau);
            R = eye(3) + 4 / (4 + theta ^ 2) * (skew(theta) + skew(theta) ^ 2 / 2);
            X = SO3_RotationMatrix(R);
        end
        function tau = inv_cayley(X)
            
        end
    end

    methods(Static)
        function X = Exp(tau)
            theta = norm(tau);
            u = tau / theta;
            R = eye(3) + sin(theta) * skew(u) + (1 - cos(theta)) * skew(u) ^ 2;
            X = SO3_RotationMatrix(R);
        end
    end
end

