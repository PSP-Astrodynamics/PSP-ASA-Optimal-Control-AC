classdef CasADi_MPC_freetf_delay_6DoF < matlab.System
    properties (Nontunable)
        vehicle = Vehicle(100000,30,60,deg2rad(20),800000,2000000, Name = "Default")
        max_t_step (1, 1) double = 0.5
        steps (1, 1) double = 50
        max_iter (1, 1) double = 400
        x_initial (12, 1) double
        %States: x, y, z, x_dot, y_dot, z_dot, theta1, theta2, theta3, w1,
        %w2, w3. z is up
        x_final (12, 1) double = [0, 0, 0, 0, 0, 0, pi/2, 0, 0, 0, 0, 0]
        delay_time (1, 1) double = 0.1
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
            num = 4;
        end
        function [dt1, dt2, dt3, dt4] = getOutputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
        end
        function [dt1, dt2] = getInputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
        end
        function [sz1, sz2, sz3, sz4] = getOutputSizeImpl(obj)
        	sz1 = [obj.steps,12];
            sz2 = [obj.steps,4];
            sz3 = [1, obj.steps];
            sz4 = [12, 1];
        end
        function [sz1, sz2] = getInputSizeImpl(~)
        	sz1 = [12, 1];
            sz2 = [12, 1];
        end
        function [cp1, cp2] = isInputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
        end
        function [cp1, cp2, cp3, cp4] = isOutputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
        end
        function [fz1, fz2] = isInputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
        end
        function [fz1, fz2, fz3, fz4] = isOutputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
        end
        function setupImpl(obj,~)
            % Define the optimization problem
            import casadi.*
            obj.opti = Opti();

            % Generate the array of state and control vectors
    
            %States: x, y, z, x_dot, y_dot, z_dot, theta1, theta2, theta3, w1,
            %w2, w3. z is up
            obj.x = obj.opti.variable(obj.steps, 12);  % steps x 12 matrix
            % Controls: thrust (i), thrust (j), thrust (k), n(rad), roll control
            obj.u = obj.opti.variable(obj.steps, 4);

            obj.xf = obj.opti.parameter(1, 12);
            obj.opti.set_value(obj.xf, obj.x_final);
            obj.opti.subject_to(obj.x(obj.steps, :) == obj.xf); % Final state


            obj.tf = obj.opti.variable(1);
                                                                  
            %Equal weighting of all control inputs, 2 * angular velocity
            %against all 3 rotational axes
            cost = sum(obj.u(:,1)) * obj.tf / obj.steps*0.2 + sum(obj.u(:,1).^2) + sum(obj.u(:,2).^2) + sum(obj.u(:,3).^2) ...
                + sum(obj.u(:,4).^2) + 2 * (sum(obj.x(:,10).^2) + sum(obj.x(:,11).^2) + sum(obj.x(:,12).^2));
            
            %cost = sum(obj.u(:,1)) * obj.tf / obj.steps * 0.2 + sum(obj.u(:,1).^2) + sum(obj.u(:,2).^2) + sum(obj.x(:,6).^2) + log(sum((obj.x(:, 1) - obj.xf(:, 1)).^2) + sum((obj.x(:, 2) - obj.xf(:, 2)).^2));
            %cost = sum(obj.u(:,1)) * obj.tf / obj.steps*0.2 + sum(obj.u(:,1))^2 * 0 + sum(obj.u(:,2).^2)*0.1 + 0*sum(obj.x(:,6).^2);
            obj.opti.minimize(cost);

            %constraints
            thrust_max = 0.10;
            gimbal_max = pi/90;
            for i = 1:(obj.steps-1)
                % Current state
                x_current = obj.x(i, :)';
                u_current = obj.u(i, :)';
                
                % Impose thrust limits
                                                                                                                            
                 total_thrust_next = sqrt((obj.u(i+1, 1).^2 + obj.u(i+1, 2).^2 + obj.u(i+1, 3).^2));
                 total_thrust = sqrt((obj.u(i, 1).^2 + obj.u(i, 2).^2 + obj.u(i, 3).^2));
                 obj.opti.subject_to((total_thrust_next - total_thrust)/(obj.tf/obj.steps) < thrust_max);
                 %obj.obti.subject_to((obj.u()))

                % %impose gimbal limits
                 %old: obj.opti.subject_to(((obj.u(i+1, 2) - obj.u(i, 2))/(obj.tf/obj.steps)) < gimbal_max);
                 gimbal_next = sqrt((obj.u(i+1, 1).^2 + obj.u(i+1, 2).^2 + obj.u(i+1, 3).^2));
                 gimbal_curr = sqrt((obj.u(i, 1).^2 + obj.u(i, 2).^2 + obj.u(i, 3).^2));
                 obj.opti.subject_to((total_thrust_next - total_thrust)/(obj.tf/obj.steps) < thrust_max);

                % Define the state derivatives
                x_dot = SymDynamicsEuler6DoF(x_current, u_current, obj.vehicle.m, obj.vehicle.L, obj.vehicle.I * [1; 1; 1], obj.vehicle.max_thrust);
                
                % Euler integration for dynamics constraints
                x_next = x_current + x_dot * obj.tf / obj.steps;
                
                % Impose the dynamics constraint
                obj.opti.subject_to(obj.x(i+1, :)' == x_next);
            end

            % Thrust percentage bounds
            %   old:
            %obj.opti.subject_to(obj.u(:,1) >= obj.vehicle.min_thrust / obj.vehicle.max_thrust);
            %obj.opti.subject_to(obj.u(:,1) <= 1);
            obj.opti.subject_to(total_thrust >= obj.vehicle.min_thrust / obj.vehicle.max_thrust);
            obj.opti.subject_to(total_thrust <= 1);
            
            % Thrust angle bounds
            %obj.opti.subject_to(obj.u(:,1) >= cos(obj.vehicle.max_gimbal) ./ vecnorm(obj.u, 2, 2));
            

            % State constraint
            obj.opti.subject_to(obj.x(:,2) >= 0);
            
            % Parameter constraints
            obj.opti.subject_to(0.5 < obj.tf)
            obj.opti.subject_to(obj.tf < obj.max_t_step * obj.steps)

            % Solver options
            p_opts = struct('expand',true,'error_on_fail',false,'verbose',false);
            s_opts = struct('max_iter',obj.max_iter);
            obj.opti.solver('ipopt', p_opts, s_opts);

            % Initial tf guess 
            [~, ~, tf_guess_initial] = guess_6DoF_with_tf(obj.x_initial', obj.x_final', obj.steps, obj.vehicle);
            obj.tf_opt = tf_guess_initial;
            % Setup parametrized initial condition
            obj.u_opt = zeros([obj.steps, 4]);
            % Propagate the initial state forward by the delay time
            x_pred = obj.propagateState(obj.x_initial);

            obj.p = obj.opti.parameter(1, 12);
            obj.opti.set_value(obj.p, x_pred');
            obj.opti.subject_to(obj.x(1, :) == obj.p); % Initial state

            % Initial guess 
            [x_guess, u_guess, tf_guess] = guess_6DoF_with_tf(obj.x_initial', obj.x_final', obj.steps, obj.vehicle);

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

        function [x_opt, u_opt, t_opt, x_pred] = stepImpl(obj, x_current, x_final)
            % Propagate the current state forward by the delay time
            [x_pred, steps_delay] = obj.propagateState(x_current);

            obj.opti.set_value(obj.xf, x_final);
            obj.opti.set_value(obj.p, x_current');
        
            % Initial guess 
            % Update the initial guess for control inputs
            % Shift the previous control inputs forward and append the last input
            obj.opti.set_initial(obj.u, [obj.u_opt((steps_delay + 1):end, :); repmat(obj.u_opt(end, :), steps_delay, 1)]);
            
            % Update the initial guess for states:
            % Shift the previous states forward and append the last state
            obj.opti.set_initial(obj.x, [x_pred'; obj.x_opt((steps_delay + 2):end, :); repmat(obj.x_opt(end, :), steps_delay, 1)]);
    
            obj.opti.set_initial(obj.tf, obj.tf_opt * (1 - steps_delay / obj.steps))

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

    methods (Access = private)
        function [x_pred, steps_delay] = propagateState(obj, x_current)

            % Calculate the number of steps corresponding to the delay time
            steps_delay = round(obj.delay_time / (obj.tf_opt / obj.steps));
            x_pred = x_current;

            for i = 1:steps_delay
                if i <= size(obj.u_opt, 1)
                    u_current = obj.u_opt(i, :)';
                else
                    u_current = obj.u_opt(end, :)';
                end

                % Calculate state derivatives
                x_dot = SymDynamicsEuler6DoF(x_pred,u_current, obj.vehicle.m, obj.vehicle.L, obj.vehicle.I * [1; 1; 1], obj.vehicle.max_thrust);

                x_pred = x_pred + x_dot * (obj.tf_opt / obj.steps);
            end
        end
    end
end