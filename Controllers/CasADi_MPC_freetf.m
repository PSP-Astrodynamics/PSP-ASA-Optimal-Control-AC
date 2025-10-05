classdef CasADi_MPC_freetf < matlab.System
    properties (Nontunable)
        vehicle = Vehicle(100000,30,60,deg2rad(20),800000,2000000, Name = "Default")
        max_t_step (1, 1) double = 0.5
        steps (1, 1) double = 50
        max_iter (1, 1) double = 400
        x_initial (6, 1) double
        x_final (6, 1) double = [0, 0, 0, 0, pi/2, 0]
    end

    properties (Access = private)
        opti
        x
        u
        tf
        xf
        p
        x_opt
        u_opt
        tf_opt
    end
    methods
        function obj = CasADi_MPC(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
    end

    methods (Access = protected)
        function num = getNumInputsImpl(~)
            num = 2;
        end
        function num = getNumOutputsImpl(~)
            num = 3;
        end
        function [dt1, dt2, dt3] = getOutputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
        end
        function [dt1, dt2] = getInputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
        end
        function [sz1, sz2, sz3] = getOutputSizeImpl(obj)
        	sz1 = [obj.steps,6];
            sz2 = [obj.steps,2];
            sz3 = [1, obj.steps];
        end
        function [sz1, sz2] = getInputSizeImpl(~)
        	sz1 = [6, 1];
            sz2 = [6, 1];
        end
        function [cp1, cp2] = isInputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
        end
        function [cp1, cp2, cp3] = isOutputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
            cp3 = false;
        end
        function [fz1, fz2] = isInputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
        end
        function [fz1, fz2, fz3] = isOutputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
            fz3 = true;
        end
        function setupImpl(obj,~)
            % Define the optimization problem
            import casadi.*
            obj.opti = Opti();

            % Generate the array of state and control vectors
    
            % States: x, y, x_dot, y_dot, theta, theta_dot
            obj.x = obj.opti.variable(obj.steps, 6);  % steps x 6 matrix
            % Controls: thrust (percent), thrust_angle (rad)
            obj.u = obj.opti.variable(obj.steps, 2);

            obj.xf = obj.opti.parameter(1, 6);
            obj.opti.set_value(obj.xf, obj.x_final);
            obj.opti.subject_to(obj.x(obj.steps, :) == obj.xf); % Final state


            obj.tf = obj.opti.variable(1);

            cost = sum(obj.u(:,1)) * obj.tf / obj.steps*0.5 + sum(obj.u(:,1).^2) + sum(obj.u(:,2).^2) + sum(obj.x(:,6).^2);
            %cost = 0*(sum(obj.u(:,1)) * obj.tf / obj.steps*0.2 + sum(obj.u(:,1))^2 * 0.5 + sum(obj.u(:,2).^2) + sum(obj.x(:,6).^2) + sum((obj.x(:, 1) - obj.xf(:, 1)).^2)) + sum((obj.x(:, 2) - obj.xf(:, 2)).^2) + 0*sum((obj.x(:, 4) - obj.xf(:, 4)).^2);
            obj.opti.minimize(cost);

            %constraints
            for i = 1:(obj.steps-1)
                % Current state
                x_current = obj.x(i, :)';
                u_current = obj.u(i, :)';
                
                % Define the state derivatives
                x_dot = Dynamics3DoF(x_current, u_current, obj.vehicle);
                
                % Euler integration for dynamics constraints
                x_next = x_current + x_dot * obj.tf / obj.steps;
                
                % Impose the dynamics constraint
                obj.opti.subject_to(obj.x(i+1, :)' == x_next);
            end

            % Thrust percentage bounds
            obj.opti.subject_to(obj.u(:,1) >= obj.vehicle.min_thrust / obj.vehicle.max_thrust);
            obj.opti.subject_to(obj.u(:,1) <= 1);
            
            % Thrust angle bounds
            obj.opti.subject_to(obj.u(:,2) >= -obj.vehicle.max_gimbal);
            obj.opti.subject_to(obj.u(:,2) <= obj.vehicle.max_gimbal);

            % State constraint
            obj.opti.subject_to(obj.x(:,2) >= 0);
            
            % Parameter constraints
            obj.opti.subject_to(0.5 < obj.tf)
            obj.opti.subject_to(obj.tf < obj.max_t_step * obj.steps)

            % Solver options
            p_opts = struct('expand',true,'error_on_fail',false,'verbose',false);
            s_opts = struct('max_iter',obj.max_iter);
            obj.opti.solver('ipopt', p_opts, s_opts);

            % Setup parametrized initial condition
            obj.p = obj.opti.parameter(1, 6);
            obj.opti.set_value(obj.p, [obj.x_initial(1), obj.x_initial(2), obj.x_initial(3), obj.x_initial(4), obj.x_initial(5), obj.x_initial(6)]);
            obj.opti.subject_to(obj.x(1, :) == obj.p); % Initial state

            % Initial guess 
            [x_guess, u_guess, tf_guess] = guess_3DoF_with_tf(obj.x_initial', obj.x_final', obj.steps, obj.vehicle);

            obj.opti.set_initial(obj.x, x_guess);
            obj.opti.set_initial(obj.u, u_guess);
            obj.opti.set_initial(obj.tf, tf_guess);

            % Solve the optimization problem
            sol = obj.opti.solve();
    
            obj.u_opt = sol.value(obj.u);
            obj.x_opt = sol.value(obj.x);
            obj.tf_opt = sol.value(obj.tf);

            % Initializing Dual Variables
            lam_g0 = sol.value(obj.opti.lam_g);
            obj.opti.set_initial(obj.opti.lam_g, lam_g0);
        end

        function [x_opt, u_opt, t_opt] = stepImpl(obj, x_current, x_final)
            obj.opti.set_value(obj.xf, x_final);
            obj.opti.set_value(obj.p, x_current'); % Should make it predict into the future to account for delay
        
            % Initial guess 
            obj.opti.set_initial(obj.u, [obj.u_opt(2:end,:); obj.u_opt(end,:)]);
            obj.opti.set_initial(obj.x, [x_current'; obj.x_opt(3:end,:); obj.x_opt(end,:)]);
            obj.opti.set_initial(obj.tf, obj.tf_opt * (1 - 1 / obj.steps))

            % Solve the optimization problem
            sol = obj.opti.solve();
    
            obj.u_opt = sol.value(obj.u);
            obj.x_opt = sol.value(obj.x);
            obj.tf_opt = sol.value(obj.tf);

            x_opt = obj.x_opt;
            u_opt = obj.u_opt;
            tf_opt = obj.tf_opt;
            t_opt = linspace(0, tf_opt, obj.steps);
        end

    end
end