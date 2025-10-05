classdef Group
    %GROUP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract)
        Name
        size
        dim
        element
        identity
    end

    methods(Abstract)
        % Abstract methods that every subclass must implement
        %evaluate(obj, heights) %- implement?

        obj = compose(X, Y)
        obj = inv(X)

        obj = constraint(X)
        obj = act(obj, x)

        tau = vee(obj, tau_hat)
        tau_hat = hat(obj, tau)

        X = Exp(obj, tau)
        tau = Log(X)

        obj = cayley(obj, X)
        obj = inv_cayley(obj, X)

        obj = Ad(X)
        obj = ad(x)
    end

    methods
        % function X = Exp(G, tau)
        %     tau_hat = G.hat(tau);
        %     X = G.exp(tau_hat);
        % end
        % function tau = Log(G, X)
        %     tau_hat = G.log(X);
        %     tau = G.vee(tau_hat);
        % end

        % Right is local frame (depends on convention)
        function Y = rplus(X, tau_X)
            Y = X.compose(X.Exp(tau_X));
        end
        function tau_X = rminus(X, Y)
            tau_X = X.inv.compose(Y).Log();
        end

        % Left is global frame (depends on convention)
        function Y = lplus(tau_eps, X)
            Y = X.Exp(tau_eps).compose(X);
        end
        function tau_eps = lminus(X, Y)
            tau_eps = Y.compose(X.inv).Log;
        end

        function [xi_r] = right_invariant_error(X_true, X_est)
            xi_r = X_est.compose(X_true.inv).Log();
        end
        function [xi_l] = left_invariant_error(X_true, X_est)
            xi_l = X_true.inv.compose(X_est).Log();
        end

        function X = sample_right_gaussian(X_mean, Sigma)
            % Sample gaussian with uncertainty in body frame
            L = chol(Sigma);
            tau = L * randn([X_mean.dim, 1]);
            X = rplus(X_mean, tau);
        end
        function X = sample_left_gaussian(Sigma, X_mean)
            % Sample gaussian with uncertainty in global frame
            L = chol(Sigma);
            tau = L * randn([X_mean.dim, 1]);
            X = lplus(tau, X_mean);
        end

        function coadjoint_map = coadj(x)
            coadjoint_map = ad(x)';
        end
    end
end

