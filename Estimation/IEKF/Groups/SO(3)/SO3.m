classdef SO3 < Group
    %SO3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim = 3
    end
    
    methods
        function obj = SO3()
            %SO3 Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

    methods(Static)
        function Adjoint = Ad(X)
            % Manifold adjoint
            Adjoint = X.R;
        end
        function adjoint = ad(G, twist)
            % Lie algebra adjoint (Lie Bracket)
            % The derivative of Ad at identity
            % Always G.hat()?
            adjoint = G.hat(twist);
        end

        function val = J_r(phi)
            % SO(3) right jacobian
            phi_norm = norm(phi);
            val = eye(3) - (1 - cos(phi_norm)) / phi_norm ^ 2 * skew(phi) + (phi_norm - sin(phi_norm)) / phi_norm ^ 3 * skew(phi) ^ 2;
        end

        function val = J_l(phi)
            % SO(3) left jacobian - transpose of right jacobian
            phi_norm = norm(phi);
            val = eye(3) + (1 - cos(phi_norm)) / phi_norm ^ 2 * skew(phi) + (phi_norm - sin(phi_norm)) / phi_norm ^ 3 * skew(phi) ^ 2;
        end

        function val = Gamma_2(phi)
            % name?
            phi_norm = norm(phi);
            val = 1/2 * eye(3) + (phi_norm - sin(phi_norm) / phi_norm ^ 3) * skew(phi) + (phi_norm ^ 2 + 2 * cos(phi_norm) - 2) / (2 * phi_norm ^ 4) * skew(phi) ^ 2;
        end

        function val = J_r_inv(theta)
            % SO(3) inverse right jacobian
            theta_norm = norm(theta);
            val = eye(3) + 1/2 * skew(theta) + (1 / theta_norm ^ 2 - (1 + cos(theta_norm)) / (2 * theta_norm * sin(theta_norm))) * skew(theta) ^ 2;
        end

        function val = J_l_inv(theta)
            % SO(3) inverse left jacobian - transpose of inverse right jacobian
            theta_norm = norm(theta);
            val = eye(3) - 1/2 * skew(theta) + (1 / theta_norm ^ 2 - (1 + cos(theta_norm)) / (2 * theta_norm * sin(theta_norm))) * skew(theta) ^ 2;
        end

        function X_dot = f_u(X, w)
            X_dot = X.compose(X.hat(w)).element;
        end 

        function X_kp1 = approx_dynamics(X_k, w, delta_t)
            X_kp1 = X_k.rplus(w * delta_t);
        end
    end
end

