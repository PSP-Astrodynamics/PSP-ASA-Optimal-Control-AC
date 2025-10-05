classdef SO3_quaternion < SO3
    %SO3_QUATERNION Summary of this class goes here
    %   Detailed explanation goes here
    %   Scalar last quaternions

      properties
        Name = "SO(3) Quaternion"
        size = 4
        element
        identity = [0; 0; 0; 1]
    end
    
    methods
        function obj = SO3_quaternion(q)
            %SO3_QUATERNION Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = q;
        end

        function val = w(obj)
            val = obj.element(4);
        end

        function val = v(obj)
            val = obj.element(1:3);
        end

        function val = R(obj)
            w = obj.w;
            v = obj.v;
            val = (w ^ 2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * skew(v);
        end
        
        function composition = compose(X, Y)
            new_quat = [X.w * Y.v + Y.w * X.v + cross(X.v, Y.v); X.w * Y.w - dot(X.v, Y.v)];
            composition = SO3_quaternion(new_quat);
        end
        function inverse = inv(X)
            inverse = SO3_quaternion([-X.v; X.w]);
        end

        function val = constraint(X)
            val = X.compose(X.inv);
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
            tau_hat = [0; theta / 2];
        end

        function X = Exp(G, tau)
            theta = norm(tau);
            u = tau / theta;
            q = [u * sin(theta / 2); cos(theta / 2)];
            X = SO3_quaternion(q);
        end
        function tau = Log(X)
            w = X.w * sign(X.w);
            v = X.v * sign(X.w);
            tau = 2 * v * atan2(norm(v), w) / norm(v);
        end

        function X = cayley(G, tau)
            theta = norm(tau);
            %R = eye(3) + 4 / (4 + theta ^ 2) * (skew(theta) + skew(theta) ^ 2 / 2);
            X = SO3_quaternion(R);
        end
        function tau = inv_cayley(X)
            
        end
    end
end

