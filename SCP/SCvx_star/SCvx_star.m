function [scvxstar_sol] = SCvx_star(prob, scvxstar_ops, parser)
%SCVX_STAR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.
% Uses augmented lagrange framework
arguments
    prob 
    scvxstar_ops 
    parser string {mustBeMember(parser, ["CVX", "CVXPY"])} = "CVX"
end

x_ref = zeros([prob.n.x, prob.N, scvxstar_ops.iter_max + 1]);
u_ref = zeros([prob.n.u, prob.Nu, scvxstar_ops.iter_max + 1]);
p_ref = zeros([prob.n.p, scvxstar_ops.iter_max + 1]);

x_ref(:, :, 1) = prob.scale_x(prob.guess.x);
u_ref(:, :, 1) = prob.scale_u(prob.guess.u);
p_ref(:, 1) = prob.scale_p(prob.guess.p);

scvxstar_sol.converged = false;
scvxstar_sol.objective = zeros([1, scvxstar_ops.iter_max + 1]);
scvxstar_sol.Delta = zeros([prob.n.x, prob.N + 1, scvxstar_ops.iter_max + 1]);
scvxstar_sol.delta_xp = zeros([1, scvxstar_ops.iter_max]);

lambda = zeros((prob.N + 1) * prob.n.x, scvxstar_ops.iter_max + 1);
mu = zeros([prob.n.ncvx_ineq_total, scvxstar_ops.iter_max + 1]);
delta = zeros([1, scvxstar_ops.iter_max + 1]);
delta(1) = inf;
w = zeros([1, scvxstar_ops.iter_max + 1]);
w(1) = scvxstar_ops.w_0;
r = zeros([1, scvxstar_ops.iter_max + 1]);
r(1) = scvxstar_ops.r_0;
rho = zeros([1, scvxstar_ops.iter_max]);
acceptance = zeros([1, scvxstar_ops.iter_max]);
delta_J = zeros([1, scvxstar_ops.iter_max]);
delta_L = zeros([1, scvxstar_ops.iter_max]);
chi = zeros([1, scvxstar_ops.iter_max]);
g_eq = zeros(size(lambda));
g_ineq = zeros(size(mu));
g_eq_relax = zeros(size(lambda));
g_ineq_relax = zeros(size(mu));

% Convexify along initial guess
[prob, scvxstar_sol.Delta(:, :, 1)] = convexify_along_reference(prob, prob.guess.x, prob.guess.u, prob.guess.p);
g_eq(1:(prob.N + 1) * prob.n.x, 1) = -reshape(scvxstar_sol.Delta(:, :, 1), [], 1);
g_eq_relax(1:(prob.N + 1) * prob.n.x, 1) = reshape(scvxstar_sol.Delta(:, :, 1), [], 1); % Not really meaningful because there is no previous linearization (disc propagation should match continuous)
[g_eq(((prob.N + 1) * prob.n.x + 1):end, 1), g_ineq(:, 1), g_eq_relax(((prob.N + 1) * prob.n.x + 1):end, 1), g_ineq_relax(:, 1)] = eval_nonconvex_constraints(prob, prob.guess.x, prob.guess.u, prob.guess.p, prob.guess.x, prob.guess.u, prob.guess.p);

disp(" k |       status      |   vd  |   vs  |  vbc_0 |  vbc_N |    J    |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   |   dp   | delta |  dyn  |  eta  | eta_x | eta_u | eta_p")

for i = 1:(scvxstar_ops.iter_max)
    % Solve convex subproblem and update reference
    [x_ref(:, :, i + 1), u_ref(:, :, i + 1), p_ref(:, i + 1), sol_info] = solve_SCvx_star_convex_subproblem(prob, scvxstar_ops, x_ref(:, :, i), u_ref(:, :, i), p_ref(:, i), lambda(:, i), mu(:, i), @(z, variant) l1_approx_penalty(z, scvxstar_ops.tau, variant), r(i), w(i));

    % Convexify along reference trajectory
    [prob, scvxstar_sol.Delta(:, :, i + 1)] = convexify_along_reference(prob, prob.unscale_x(x_ref(:, :, i + 1)), prob.unscale_u(u_ref(:, :, i + 1)), prob.unscale_p(p_ref(:, i + 1)));

    % Compute original and convexified constraint violations
    g_eq(1:(prob.N + 1) * prob.n.x, i + 1) = -reshape(scvxstar_sol.Delta(:, :, i + 1), [], 1);
    [g_eq(((prob.N + 1) * prob.n.x + 1):end, i + 1), g_ineq(:, i + 1), t1, t2] = eval_nonconvex_constraints(prob, x_ref(:, :, i + 1), u_ref(:, :, i + 1), p_ref(:, i + 1), x_ref(:, :, i), u_ref(:, :, i), p_ref(:, i));
    g_eq_relax(:, i + 1) = sol_info.xi;
    g_ineq_relax(:, i + 1) = sol_info.zeta;

    % Compute original and convexified objectives for solution
    J_NL_aug = prob.objective(prob.unscale_x(x_ref(:, :, i + 1)), prob.unscale_u(u_ref(:, :, i + 1)), prob.unscale_p(p_ref(:, i + 1))) + penalty_function(lambda(:, i), g_eq(:, i + 1), mu(:, i), g_ineq(:, i + 1), @(z, variant) l1_approx_penalty(z, scvxstar_ops.tau, variant), w(i));
    J_NL_aug_ref = prob.objective(prob.unscale_x(x_ref(:, :, i)), prob.unscale_u(u_ref(:, :, i)), prob.unscale_p(p_ref(:, i))) + penalty_function(lambda(:, i), g_eq(:, i), mu(:, i), g_ineq(:, i), @(z, variant) l1_approx_penalty(z, scvxstar_ops.tau, variant), w(i));
    J_aug = sol_info.J + sol_info.J_vc;

    % Update algorithm weights
    [rho(i), acceptance(i), delta_J(i), delta_L(i)] = solution_acceptance_criterion(J_NL_aug, J_NL_aug_ref, J_aug, scvxstar_ops.eta_0);

    rho(i)

    if acceptance(i)
        if abs(delta_J(i)) < delta(i)
            [lambda(:, i + 1), mu(:, i + 1), w(:, i + 1)] = multiplier_update(lambda(:, i), mu(:, i), w(:, i), g_eq_relax(:, i + 1), g_ineq_relax(:, i + 1), scvxstar_ops.tau, scvxstar_ops.beta, scvxstar_ops.w_max);
            delta(i + 1) = stationary_tolerance_update(delta(i), delta_J(i), scvxstar_ops.gamma);
        end
    else
        % Don't accept solution update
        x_ref(:, :, i + 1) = x_ref(:, :, i);
        u_ref(:, :, i + 1) = u_ref(:, :, i);
        p_ref(:, i + 1) = p_ref(:, i);
    end

    if ~acceptance(i) || abs(delta_J(i)) >= delta(i)
        lambda(:, i + 1) = lambda(:, i);
        mu(:, i + 1) = mu(:, i);
        w(i + 1) = w(i);
        delta(i + 1) = delta(i);
    end

    r(i + 1) = trust_region_update(rho(i), r(i), scvxstar_ops.eta_1, scvxstar_ops.eta_2, scvxstar_ops.alpha_1, scvxstar_ops.alpha_2, scvxstar_ops.r_min, scvxstar_ops.r_max);

    % Check stopping criteria (4.30) - FOR COMPARISON AGAINST PTR PURPOSES
    scvxstar_sol.delta_xp(i) = ptr_stopping(x_ref(:, :, i + 1), p_ref(:, i + 1), x_ref(:, :, i), p_ref(:, i), 2);

    % Display results of iteration
    sol_info.delta = scvxstar_sol.delta_xp(i);
    sol_info.dyn = sum(vecnorm(scvxstar_sol.Delta(:, :, i + 1), 2, 1)) <= scvxstar_ops.feas_tol;
    % CHECK CONSTRAINT INFEASIBILITY TOO
    if sol_info.dyn
        dyn_str = string(sol_info.dyn);
    else
        dyn_str = sum(vecnorm(scvxstar_sol.Delta(:, :, i + 1), 2, 1));
    end
    chi(i) = norm([g_eq(:, i + 1); max(g_ineq(:, i + 1), 0)]);

    scvxstar_sol.info(i) = sol_info;

    if i == 1
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g |         | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, scvxstar_sol.info(i).status, norm(scvxstar_sol.info(i).vd), norm(scvxstar_sol.info(i).vs), norm(scvxstar_sol.info(i).vbc_0), norm(scvxstar_sol.info(i).vbc_N), scvxstar_sol.info(i).J, scvxstar_sol.info(i).J_tr, scvxstar_sol.info(i).J_vc, sum(scvxstar_sol.info(i).dx), sum(scvxstar_sol.info(i).du), sum(scvxstar_sol.info(i).dp), scvxstar_sol.info(i).delta, dyn_str, norm(scvxstar_sol.info(i).eta), norm(scvxstar_sol.info(i).eta_x), norm(scvxstar_sol.info(i).eta_u), norm(scvxstar_sol.info(i).eta_p))
    else
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g | %6.3f | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, scvxstar_sol.info(i).status, norm(scvxstar_sol.info(i).vd), norm(scvxstar_sol.info(i).vs), norm(scvxstar_sol.info(i).vbc_0), norm(scvxstar_sol.info(i).vbc_N), scvxstar_sol.info(i).J, scvxstar_sol.info(i).J_tr, scvxstar_sol.info(i).J_vc, scvxstar_sol.info(i).dJ, sum(scvxstar_sol.info(i).dx), sum(scvxstar_sol.info(i).du), sum(scvxstar_sol.info(i).dp), scvxstar_sol.info(i).delta, dyn_str, norm(scvxstar_sol.info(i).eta), norm(scvxstar_sol.info(i).eta_x), norm(scvxstar_sol.info(i).eta_u), norm(scvxstar_sol.info(i).eta_p))
    end

    % Check stopping conditions
    if i >= scvxstar_ops.iter_min && delta_J(i) <= scvxstar_ops.opt_tol && ~scvxstar_sol.converged && chi(i) <= scvxstar_ops.feas_tol %sol_info.dyn
        scvxstar_sol.converged = true;
        scvxstar_sol.converged_i = i;

        break
    end
end

if scvxstar_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, scvxstar_sol.delta_xp(i), sum(vecnorm(scvxstar_sol.Delta(:, :, end), 2, 1)))
end

scvxstar_sol.lambda = lambda;
scvxstar_sol.delta = delta;
scvxstar_sol.delta = delta;
scvxstar_sol.w = w;
scvxstar_sol.r = r;
scvxstar_sol.rho = rho;
scvxstar_sol.acceptance = acceptance;
scvxstar_sol.delta_J = delta_J;
scvxstar_sol.delta_L = delta_L;
scvxstar_sol.chi = chi;
scvxstar_sol.g_ineq = g_ineq;
scvxstar_sol.g_eq_relax = g_eq_relax;
scvxstar_sol.g_ineq_relax = g_ineq_relax;

scvxstar_sol.x = prob.unscale_x(x_ref);
scvxstar_sol.u = prob.unscale_u(u_ref);
scvxstar_sol.p = prob.unscale_p(p_ref);
end

function [rho, acceptance, delta_J, delta_L] = solution_acceptance_criterion(J_NL_aug, J_NL_aug_ref, J_aug, eta_0)
    delta_J = J_NL_aug_ref - J_NL_aug; % Nonlinear cost reduction
    delta_L = J_NL_aug_ref - J_aug; % Approximated cost reduction

    if delta_L == 0
        rho = 1;
    else
        rho = delta_J / delta_L;
    end
    
    acceptance = 1 - eta_0 <= rho && rho <= 1 + eta_0;
end

function [r_update] = trust_region_update(rho, r, eta_1, eta_2, alpha_1, alpha_2, r_min, r_max)
    if 1 - eta_2 <= rho && rho <= 1 + eta_2
        r_update = max(alpha_2 * r, r_min);
    elseif 1 - eta_1 <= rho && rho <= 1 + eta_1
        r_update = r;
    else
        r_update = min(r / alpha_1, r_max);
    end
end

function [lambda_update, mu_update, w_update] = multiplier_update(lambda, mu, w, g_eq_relax, g_ineq_relax, tau, beta, w_max)
    % MULTIPLIER_UPDATE Update lagrange multipliers for convexified
    % constraints
    lambda_update = lambda + l1_approx_penalty_gradient(w * g_eq_relax, tau);
    mu_update = max(mu + l1_approx_penalty(w * g_ineq_relax, tau), 0);
    w_update = min(beta * w, w_max); % update penalty weight
end

function [delta_update] = stationary_tolerance_update(delta, delta_J, gamma)
    if delta == inf
        delta_update = abs(delta_J);
    else
        delta_update = gamma * delta;
    end
end

function [delta_xp] = ptr_stopping(X, p, x_ref, p_ref, q)
    delta_xp = zero_if_empty(vecnorm(p - p_ref, q, 1)) + max(vecnorm(X - x_ref, q, 1));
end