function [x_sol, u_sol, p_sol, sol_info] = solve_SCvx_star_convex_subproblem(prob, scvxstar_ops, x_ref, u_ref, p_ref, lambda, mu, phi, r, w)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

t_k = linspace(0, prob.tf, prob.N);

cvx_begin quiet
    variable X(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable p(prob.n.p, 1)
    variable V(prob.n.x, prob.N - 1)
    %variable xi(prob.n.ncvx_eq)
    variable zeta(numel(mu))
    variable v_0(prob.n.x, 1)
    variable v_N(prob.n.x, 1)
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p)) ...%, prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref)) ...
        + penalty_function(lambda, [v_0; V(:); v_N], mu, zeta, phi, w) )
    subject to
        % Dynamics
        if prob.u_hold == "ZOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k));
            end
        elseif prob.u_hold == "FOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_minus_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.B_plus_k(:, :, k) * prob.unscale_u(U(:, k + 1)) ...
                             + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k));
            end
        end

        ncvx_ineq_index = 0; 
        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                cc_k = prob.convex_constraints{cc}{1};
                if ismember(k, cc_k)
                    cvx_constraint_func = prob.convex_constraints{cc}{2};
                    cvx_constraint_func(t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), prob.unscale_p(p)) <= 0;
                end
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                nc_k = prob.nonconvex_constraints{nc}{1};
                if ismember(k, nc_k)
                    ncvx_ineq_index = ncvx_ineq_index + 1;

                    ncvx_constraint_func = prob.nonconvex_constraints{nc}{2};
                    ncvx_constraint_func(t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), prob.unscale_p(p), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref), k) ...
                        - zeta(ncvx_ineq_index) <= 0;
                end
            end
        end
        zeta >= 0;

        % Boundary Conditions
        prob.initial_bc(prob.unscale_x(X(:, 1)), prob.unscale_p(p)) + v_0 == 0;
        prob.terminal_bc(prob.unscale_x(X(:, prob.N)), prob.unscale_p(p)) + v_N == 0;

        % Trust Region Constraints
        norm(scvxstar_ops.D_x * (X(:) - x_ref(:)), inf) <= r;
        norm(scvxstar_ops.D_u * (U(:) - u_ref(:)), inf) <= r;
        norm(scvxstar_ops.D_p * (p(:) - p_ref(:)), inf) <= r;
cvx_end

x_sol = X;
u_sol = U;
p_sol = p;

sol_info.status = cvx_status;
sol_info.vd = V;
sol_info.vs = zeta;
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p));%, prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref));
sol_info.J_tr = 0;
sol_info.J_vc = penalty_function(lambda, [V(:); v_0; v_N], mu, zeta, phi, w);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p)) ..., prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref)) ...
    - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref))) ..., prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref))) ...
    / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref));..., prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref));
sol_info.dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), 2, 1);
sol_info.du = vecnorm(U - u_ref, 2, 1);
sol_info.dp = vecnorm(p - p_ref, 2, 1);
sol_info.eta = r;
sol_info.eta_x = scvxstar_ops.D_x * norm(X(:) - x_ref(:), inf);
sol_info.eta_u = scvxstar_ops.D_u * norm(U(:) - u_ref(:), inf);
sol_info.eta_p = scvxstar_ops.D_p * norm(p(:) - p_ref(:), inf);

sol_info.xi = [v_0; V(:); v_N];
sol_info.zeta = zeta;
end