classdef TFG_RotationMatrix < Group
    %SE2_3_ROTATIONMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "Two Frame Group Rotation Matrix"
        size = 15
        dim = 9
        element
        identity = eye(5)
    end
    
    methods
        function obj = SE2_3_RotationMatrix(R, v, r)
            %SE2_3_ROTATIONMATRIX Construct an instance of this class
            %   Detailed explanation goes here
            obj.element = [R, v, r; zeros([2, 3]), eye(2)];
        end

        function val = R(obj)
            val = obj.element(1:3, 1:3);
        end
        function val = r(obj)
            val = obj.element(1:3, 5);
        end
        function val = v(obj)
            val = obj.element(1:3, 4);
        end
        
        function composition = compose(X, Y)
            % composition = X.element * Y.element - just matrix multiplication
            composition = SE2_3_RotationMatrix(X.R * Y.R, X.R * Y.v + X.v, X.R * Y.r + X.r);
        end
        function inverse = inv(X)
            inverse = SE2_3_RotationMatrix(X.R', -X.R' * X.v, -X.R' * X.r);
        end

        function val = constraint(X)
            val = X.R' * X.R;
        end
        function val = act(X, p)
            val = X.r + X.R * p;
        end

        function tau = vee(G, tau_hat)
            theta = SO3_RotationMatrix.vee(tau_hat(1:3, 1:3));
            v = tau_hat(1:3, 4);
            rho = tau_hat(1:3, 5);
            tau = [theta; v; rho];
        end
        function tau_hat = hat(G, tau)
            theta = tau(1:3);
            v = tau(4:6);
            rho = tau(7:9);
            tau_hat = [skew(theta), v, rho; 
                       zeros(2, 5)];
        end

        function X = Exp(G, tau)
            theta = tau(1:3);
            tang_v = tau(4:6);
            rho = tau(7:9);
            % Component-wise way
            %R = SO3_RotationMatrix.Exp(theta).R;
            %v = SO3.J_l(theta) * tang_v;
            %r = SO3.J_l(theta) * rho;
            %X = SE2_3_RotationMatrix(R, v, r);

            % Classic way - equivalent to component-wise way
            S = G.hat(tau);
            tau_norm = norm(tau);
            X_exp = eye(5) + S + (1 - cos(tau_norm)) / tau_norm ^ 2 * S ^ 2 + (tau_norm - sin(tau_norm)) / tau_norm ^ 3 * S ^ 3;
            
            X = G;
            X.element = X_exp;
        end
        function tau = Log(X)
            SO3_R = SO3_RotationMatrix(X.R);
            theta = SO3_R.Log();
            v = SO3.J_l_inv(theta) * X.v;
            rho = SO3.J_l_inv(theta) * X.r;
            tau = [theta; v; rho];
        end

        function obj = cayley(obj, X)

        end
        function obj = inv_cayley(obj, X)

        end

        function Adjoint = Ad(X)
            % Manifold adjoint
            Adjoint = [X.R, zeros([3, 6]);
                       skew(X.v) * X.R, X.R, zeros(3,3);
                       skew(X.r) * X.R, zeros(3, 3), X.R];
        end
        function adjoint = ad(G, twist)
            % Tangent space adjoint
            %adjoint = Ad(G.Exp(twist));
            twist_theta = twist(1:3);
            twist_v = twist(4:6);
            twist_r = twist(7:9);
            adjoint = [skew(twist_theta), zeros([3, 6]);
                       skew(twist_v), skew(twist_theta), zeros(3, 3);
                       skew(twist_r), zeros(2, 2), skew(twist_theta)];
        end

        function val = Gamma_0(X, phi)
            R = SO3_RotationMatrix(X.R);
            val = R.Exp(phi).element;
        end

        function val = Gamma_1(X, phi)
            % SO(3) left jacobian
            phi_norm = norm(phi);
            val = eye(3) ;%+ (1 - cos(phi_norm)) / phi_norm ^ 2 * skew(phi) + (phi_norm - sin(phi_norm)) / phi_norm ^ 3 * skew(phi) ^ 2;
        end
        function val = Gamma_2(obj, phi)
            phi_norm = norm(phi);
            val = 1/2 * eye(3);% + (phi_norm - sin(phi_norm)) / phi_norm ^ 3 * skew(phi) + (phi_norm ^ 2 + 2 * cos(phi_norm) - 2) / (2 * phi_norm ^ 4) * skew(phi) ^ 2;
        end

        function val = V_inv(obj, theta)
            % SO(3) inverse left jacobian
            theta_norm = norm(theta);
            val = eye(3) - 1/2 * skew(theta) + (1 / theta_norm ^ 2 - (1 + cos(theta_norm)) / (2 * theta_norm * sin(theta_norm))) * skew(theta) ^ 2;
        end

        function X_dot = f_u(X, a, w, g)
            X_dot = [X.R * skew(w), X.R * a + g, X.v;
                     zeros(2, 5)];
        end

        function val = A_r(g)
            % Log Linear Right Invariant Dynamics
            val = [zeros([3, 9]); 
                   skew(g), zeros(3, 6); 
                   zeros([3, 3]), eye(3), zeros([3, 3])];
        end

        function STM = STM_r(g, dt)
            STM = eye(9) + A_r(g) * dt + [zeros([6, 9]); 1/2 * skew(g) * dt ^ 2, zeros(3, 6)]; 
        end

        function val = A_l(u)
            % Log Linear Left Invariant Dynamics
            u_w = u(1:3);
            u_a = u(4:6);
            val = [-skew(u_w), zeros([3, 6]); 
                   -skew(u_a), -skew(u_w), zeros(3, 3); 
                   zeros([3, 3]), eye(3), -skew(u_w)];
        end

        function STM = STM_l(G, u, dt)
            a = u(1:3);
            w = u(4:6);

            G_R = SO3_RotationMatrix(G.R);

            STM = zeros([9, 9]);
            STM(1:3, 1:3) = G.Gamma_0(w * dt)';
            STM(4:6, 1:3) = -STM(1:3, 1:3) * G_R.hat(G.Gamma_1(w * dt) * a) * dt; 
            STM(7:9, 1:3) = -STM(1:3, 1:3) * G_R.hat(G.Gamma_2(w * dt) * a) * dt ^ 2;
            STM(4:6, 4:6) = STM(1:3, 1:3);
            STM(7:9, 4:6) = STM(1:3, 1:3) * dt;
            STM(7:9, 7:9) = STM(1:3, 1:3);
        end

        function X_kp1 = approx_dynamics(X_k, u, g, delta_t, true)
            % Exact if constant IMU measurements over delta_t
            a = u(1:3);
            w = u(4:6);

            SO3_kp1 = SO3_RotationMatrix(X_k.R);
            R_kp1 = SO3_kp1.rplus(w * delta_t).element;
            S3_kp1 = SO3_quaternion([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0] * dcm2quat(X_k.R)');
            dR = quat_rotmatrix(S3_kp1.Exp(w * delta_t).element);
            q_kp1 = S3_kp1.rplus(w * delta_t).element;
            R_kp1 = quat_rotmatrix(q_kp1);
            v_kp1 = X_k.v + X_k.R * X_k.Gamma_1(w * delta_t) * a * delta_t + g * delta_t;
            r_kp1 = X_k.r + X_k.v * delta_t + X_k.R * X_k.Gamma_2(w * delta_t) * a * delta_t ^ 2 + 1 / 2 * g * delta_t ^ 2;

            X_kp1 = SE2_3_RotationMatrix(R_kp1, v_kp1, r_kp1);
        end
    end
end


function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);
    x = v(1);
    y = v(2);
    z = v(3);
    R = [w ^ 2 + x ^ 2 - y ^ 2 - z ^ 2, 2 * (x * y - w * z), 2 * (x * z + w * y);
         2 * (x * y + w * z), w ^ 2 - x ^ 2 + y ^ 2 - z ^ 2, 2 * (x * z - w * x);
         2 * (x * z - w * y), 2 * (y * z + w * x), w ^ 2 - x ^ 2 - y ^ 2 + z ^ 2]';
end