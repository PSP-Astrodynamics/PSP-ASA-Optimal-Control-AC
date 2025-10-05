function [x_sol, u_sol, sol_info, problem] = solve_ptr_convex_subproblem_no_p_CVXPY(prob, ptr_ops, x_ref, u_ref, problem)
%SOLVE_PTR_CONVEX_SUBPROBLEM_NO_P_CVXPY Summary of this function goes here
%   Detailed explanation goes here

t1 = tic;

% Solve
if prob.Name == "Deterministic_6DoF_fixed"
    if prob.u_hold == "ZOH"
        [X_py, U_py, eta_py, V_py, v_0_py, v_N_py, solve_status_py, problem] = pyrunfile("Deterministic_6DoF_fixed_ZOH.py", ["X_sol", "U_sol", "eta", "V", "v_0", "v_N", "solve_status", "problem"], x_ref = x_ref, u_ref = u_ref, x_0 = prob.x0, x_f = prob.xf, A_k = prob.disc.A_k, B_k = prob.disc.B_k, c_k = prob.disc.c_k, params = prob.params, N = prob.N - 1, delta_t = prob.tf / (prob.N - 1), w_vc = ptr_ops.w_vc, w_tr = ptr_ops.w_tr, problem = problem);
    elseif prob.u_hold == "FOH"
        [X_py, U_py, eta_py, V_py, v_0_py, v_N_py, solve_status_py, problem] = pyrunfile("Deterministic_6DoF_fixed_FOH.py", ["X_sol", "U_sol", "eta", "V", "v_0", "v_N", "solve_status", "problem"], x_ref = x_ref, u_ref = u_ref, x_0 = prob.x0, x_f = prob.xf, A_k = prob.disc.A_k, B_k_minus = prob.disc.B_minus_k, B_k_plus = prob.disc.B_plus_k, c_k = prob.disc.c_k, params = prob.params, N = prob.N - 1, delta_t = prob.tf / (prob.N - 1), w_vc = ptr_ops.w_vc, w_tr = ptr_ops.w_tr, problem = problem);
    end
end

% Convert solution to Matlab friendly types
X = double(X_py);
U = double(U_py);
eta = double(eta_py);
V = double(V_py);
v_prime = [];
v_0 = double(v_0_py);
v_N = double(v_N_py);

solve_status = string(solve_status_py);

t2 = toc(t1);

% Package outputs

x_sol = X;
u_sol = U;
%fprintf("CVXPY Time: %.3f ms\n", t2 * 1000 )
sol_info.status = solve_status;
sol_info.vd = V;
sol_info.vs = v_prime;
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0);
sol_info.J_tr = trust_region_cost(eta, 0, ptr_ops.w_tr, 0);
sol_info.J_vc = virtual_control_cost(V, v_prime, v_0, v_N, ptr_ops.w_vc);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0) - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0)) / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0);
sol_info.dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1);
sol_info.du = vecnorm(U - u_ref, ptr_ops.q, 1);
sol_info.dp = 0;
sol_info.eta = eta;
sol_info.eta_x = 0;
sol_info.eta_u = 0;
sol_info.eta_p = 0;
end

function [J_tr] = trust_region_cost(eta, eta_p, w_tr, w_tr_p)
    J_tr = w_tr * eta' + w_tr_p * eta_p;
end

function [J_vc] = virtual_control_cost(V, v_prime, v_0, v_N, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + sum(norms(V, 1, 1)) + norm(v_0, 1) + norm(v_N, 1));
end

