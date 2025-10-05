classdef SE2_2_RotationMatrix < Group
    %SE2_2_ROTATIONMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "SE2(2) Rotation Matrix"
        size = 16
        dim = 5
        element
        identity = eye(5)
    end
    
    methods
        function obj = SE2_2_RotationMatrix(R, r, v)
            %SE2_2_ROTATIONMATRIX Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = [R, v, r; zeros([2, 1]), eye(2)];
        end

        function val = R(obj)
            val = obj.element(1:2, 1:2);
        end
        function val = r(obj)
            val = obj.element(1:2, 4);
        end
        function val = v(obj)
            val = obj.element(1:2, 3);
        end
        
        function composition = compose(X, Y)
            % composition = X.element * Y.element - just matrix multiplication
            composition = SE2_2_RotationMatrix(X.R * Y.R, X.R * Y.v + X.v, X.R * Y.r + X.r);
        end
        function inverse = inv(X)
            inverse = [X.R', -X.r' * X.v, -X.R' * X.r;
                       zeros(2, 1), eye(2)];
        end

        function val = constraint(X)
            val = X.R' * X.R;
        end
        function val = act(X, p)
            val = X.r + X.R * p;
        end

        function tau = vee(tau_hat)
            theta = tau_hat(2, 1);
            v = tau_hat(1:2, 3);
            rho = tau_hat(1:2, 4);
            tau = [theta; v; rho];
        end
        function tau_hat = hat(tau)
            theta = tau(1);
            v = tau(2:3);
            rho = tau(4:5);
            tau_hat = [skew2(theta), v, rho; 
                       zeros(2, 3)];
        end

        function X = Exp(G, tau)
            theta = tau(1);
            tang_v = tau(2:3);
            rho = tau(4:5);
            R = skew(theta);
            v = tang_v;
            r = G.V(theta) * rho;
            X = SE2_RotationMatrix(R, v, r);
        end
        function tau = Log(X)
            theta = arctan2(X.R(2, 1), X.R(1, 1));
            v = X.v;
            rho = V_inv(theta) * X.r;
            tau = [theta; v; rho];
        end

        function obj = cayley(obj, X)

        end
        function obj = inv_cayley(obj, X)

        end

        function Adjoint = Ad(X)
            % Manifold adjoint
            Adjoint = [1, zeros([1, 4]);
                       -star(X.v), X.R, zeros(2,2);
                       -star(X.r), zeros(2, 2); X.R];
        end
        function adjoint = ad(G, twist)
            % Tangent space adjoint
            %adjoint = Ad(G.Exp(twist));
            twist_theta = twist(1);
            twist_v = twist(2:3);
            twist_r = twist(4:5);
            adjoint = [zeros([1, 5]);
                       -star(twist_v), skew2(twist_theta), zeros(2, 2);
                       -star(twist_r), zeros(2, 2), skew2(twist_theta)];
        end

        function val = V(obj, theta)
            % SO(2) left jacobian
            val = sin(theta) / theta * eye(2) + (1 - cos(theta)) / theta * skew2(1);
        end
        function val = V_inv(obj, theta)
            % SO(2) inverse left jacobian
            val = sin(theta) / theta * eye(2) + (1 - cos(theta)) / theta * skew2(1);
        end

        function xi_star = star(xi)
            xi_star = [-xi(2); xi(1)];
        end

        function X_dot = f_u(X, a, w, g)
            X_dot = [X.R * skew2(w), g + X.R * a, X.v;
                     zeros(2, 4)];
        end

        function X_kp1 = approx_dynamics(X_k, u, w, g, delta_t)
            R_kp1 = X_k.R * expm(delta_t * skew2(w));
            v_kp1 = X_k.v + delta_t * (g + X_k.R * u);
            r_kp1 = X_k.r + delta_t * X_k.v;

            X_kp1 = SE2_2_RotationMatrix(R_kp1, r_kp1, v_kp1);
        end

        function A_t = A(g)
            A_t = [zeros(2, 2), zeros(2, 2), zeros(2, 2);
                   skew2(g), zeros(2, 2), zeros(2, 2);
                   zeros(2, 2), eye(2), zeros(2, 2)];
        end
    end
end

