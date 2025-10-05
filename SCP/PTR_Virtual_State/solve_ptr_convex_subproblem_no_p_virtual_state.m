function [x_sol, u_sol, sol_info] = solve_ptr_convex_subproblem_no_p_virtual_state(prob, ptr_ops, x_ref, u_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

%t1 = tic;
t_k = linspace(0, prob.tf, prob.N);

cvx_begin quiet
    variable X(prob.n.x, prob.N)
    variable xi(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable v_prime(prob.n.ncvx)
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0) ...
        + virtual_state_cost(X, xi, ptr_ops.w_vse) ...
        + trust_region_cost(X, U, 0, x_ref, u_ref, 0, ptr_ops.w_tr, 0) ...
        + sum(ptr_ops.w_prime .* v_prime))
    subject to
        % Dynamics
        if prob.u_hold == "ZOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.c_k(:, :, k));
            end
        elseif prob.u_hold == "FOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_minus_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.B_plus_k(:, :, k) * prob.unscale_u(U(:, k + 1)) ...
                             + prob.disc.c_k(:, :, k));
            end
        end

        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                cc_k = prob.convex_constraints{cc}{1};
                if ismember(k, cc_k)
                    cvx_constraint_func = prob.convex_constraints{cc}{2};
                    cvx_constraint_func(t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0) <= 0;
                end
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                nc_k = prob.nonconvex_constraints{nc}{1};
                if ismember(k, nc_k)
                    ncvx_constraint_func = prob.nonconvex_constraints{nc}{2};
                    ncvx_constraint_func(t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0, k) ...
                        - v_prime(nc) <= 0;
                end
            end
        end
        v_prime >= 0;

        % Boundary Conditions
        %prob.initial_bc(prob.unscale_x(X(:, 1)), 0) == 0;
        prob.initial_bc(prob.unscale_x(xi(:, 1)), 0) == 0;
        prob.terminal_bc(prob.unscale_x(xi(:, prob.N)), 0, prob.unscale_x(prob.unscale_x(x_ref(:, prob.N))), 0) == 0;
cvx_end

%t2 = toc(t1);

if size(U,1) == 5
    figure
    
    tiledlayout(1, 3)
    nexttile
    plot3(X(1, :), X(2, :), X(3, :)); 
    grid on
    axis equal
    
    nexttile
    stairs(U(1, :)); hold on
    stairs(U(2, :)); hold on
    stairs(U(3, :)); hold on
    stairs(U(4, :)); hold off
    grid on
    
    nexttile
    stairs(rad2deg(U(5, :))); hold on
    stairs(acosd(U(1, :) ./ U(4, :)));
    grid on
end

x_sol = X;
u_sol = U;
%fprintf("Mosek Time: %.3f ms", t2 * 1000 )
sol_info.status = cvx_status;
sol_info.vd = X - xi;
sol_info.vs = v_prime;
sol_info.vbc_0 = norm(X(:, 1) - xi(:, 1));
sol_info.vbc_N = norm(X(:, end) - xi(:, end));
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0);
sol_info.J_tr = trust_region_cost(X, U, 0, x_ref, u_ref, 0, ptr_ops.w_tr, 0);
sol_info.J_vse = virtual_state_cost(X, xi, ptr_ops.w_vse);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0) - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0)) / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0);
sol_info.dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1);
sol_info.du = vecnorm(U - u_ref, ptr_ops.q, 1);
sol_info.dp = 0;
sol_info.eta = vecnorm(X - xi, 2, 1);
sol_info.eta_x = 0;
sol_info.eta_u = 0;
sol_info.eta_p = 0;
end

% function [J_tr] = trust_region_cost(x, u, p, x_ref, u_ref, p_ref, w_tr, w_tr_p)
%     J_tr = sum(w_tr * (norms(x - x_ref, 2, 1) + norms(u - u_ref, 2, 1))) + w_tr_p * norms(p - p_ref, 2, 1);
% end
% 
% function [J_vse] = virtual_state_cost(x, xi, w_vse)
%     J_vse = sum(w_vse * norms(x - xi, 2, 1));
% end

% 2-norm squared versions (seem to result in 2 times more ptr iterations than just 2-norm)
function [J_tr] = trust_region_cost(x, u, p, x_ref, u_ref, p_ref, w_tr, w_tr_p)
    J_tr = w_tr * (sum(sum_square(x - x_ref)) + sum(sum_square(u - u_ref))) + w_tr_p * sum_square(p - p_ref);
end

function [J_vse] = virtual_state_cost(x, xi, w_vse)
    J_vse = sum(w_vse * sum_square(x - xi));
end