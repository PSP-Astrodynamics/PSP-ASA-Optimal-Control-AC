                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                classdef CasADi_MPC_nDoF < matlab.System
    properties (Nontunable)
        vehicle = Vehicle(100000,30,60,deg2rad(20),800000,2000000, Name = "Default")

        dof (1, 1)% {mustBeA(dof, "DoF")}

        t_step (1, 1) double = 0.01
        steps (1, 1) double = 400
        max_iter (1, 1) double = 100
        x_initial (:, 1) double
        x_final (:, 1) double

        cost_x (:, 1) double
        cost_u (:, 1) double
    end

    properties (Access = private)
        opti
        x
        u
        xf
        p
        x_opt
        u_opt
    end
    methods
        function obj = CasADi_MPC_nDoF(varargin)
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
        	sz1 = [obj.steps, obj.dof.nx];
            sz2 = [obj.steps, obj.dof.nu];
            sz3 = [1, obj.steps];
        end
        function [sz1, sz2] = getInputSizeImpl(~)
        	sz1 = [obj.dof.nx, 1];
            sz2 = [obj.dof.nx, 1];
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
            obj.x = obj.opti.variable(obj.steps, obj.dof.nx);  % steps x 6 matrix
            % Controls: thrust (percent), thrust_angle (rad)
            obj.u = obj.opti.variable(obj.steps, obj.dof.nu);

            obj.xf = obj.opti.parameter(1, obj.dof.nx);
            obj.opti.set_value(obj.xf, obj.x_final);
            obj.opti.subject_to(obj.x(obj.steps, 3) == obj.xf(3)); % Final state

            path_cost_x = sum(obj.cost_x' .* sum(obj.x.^2, 1));
            path_cost_u = sum(obj.cost_u' .* sum(obj.u.^2, 1));
            cost = path_cost_x + path_cost_u;

            obj.opti.minimize(cost);

            %constraints
            for i = 1:(obj.steps-1)
                % Current state
                x_current = obj.x(i, :)';
                u_current = obj.u(i, :)';
                
                % Define the state derivatives
                x_dot = Dynamics_nDoF(x_current, u_current, obj.vehicle, obj.dof);
                
                % Euler integration for dynamics constraints
                x_next = x_current + x_dot * obj.t_step;
                
                % Impose the dynamics constraint
                obj.opti.subject_to(obj.x(i+1, :)' == x_next);
            end

            % Thrust percentage bounds
            obj.opti.subject_to(sum(obj.u(:,obj.dof.ithrust).^2, 2) >= (obj.vehicle.min_thrust / obj.vehicle.max_thrust).^2);
            obj.opti.subject_to(sum(obj.u(:,obj.dof.ithrust).^2, 2) <= 1);

            % Add other constraints including thrust angle bounds
            for c = 1:numel(obj.dof.constraints)
                constraint = obj.dof.constraints{c};
                obj.opti.subject_to(constraint(obj.x, obj.u, obj.vehicle));
            end

            % State constraint
            obj.opti.subject_to(obj.x(:,obj.dof.ir(end)) >= 0); % - replace with glideslope?
            
            % Solver options
            p_opts = struct('expand',true,'error_on_fail',false,'verbose',false);
            s_opts = struct('max_iter',obj.max_iter);
            obj.opti.solver('ipopt', p_opts, s_opts);

            % Setup parametrized initial condition
            obj.p = obj.opti.parameter(1, obj.dof.nx);
            obj.opti.set_value(obj.p, obj.x_initial);
            obj.opti.subject_to(obj.x(1, :) == obj.p); % Initial state

            % Initial guess 
            [x_guess, u_guess] = guess_nDoF(obj.x_initial', obj.x_final', obj.steps, obj.t_step, obj.vehicle, obj.dof);

            obj.opti.set_initial(obj.x, x_guess);
            obj.opti.set_initial(obj.u, u_guess);

            % Solve the optimization problem
            sol = obj.opti.solve();
    
            obj.u_opt = sol.value(obj.u);
            obj.x_opt = sol.value(obj.x);

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

            % Solve the optimization problem
            sol = obj.opti.solve();
    
            u_opt = sol.value(obj.u);
            x_opt = sol.value(obj.x);
            t_opt = linspace(0, obj.t_step * obj.steps, obj.steps);

            obj.u_opt = u_opt;
            obj.x_opt = x_opt;
        end

    end
end