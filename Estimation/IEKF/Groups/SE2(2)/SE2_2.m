classdef SE2_2 < Group
    %SE2_2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "SE2(2)"
        size
        dim
        element
        identity
    end
    
    methods
        function obj = SE2_2(inputArg1,inputArg2)
            %SE2_2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function  = compose(obj, X, Y)
            
        end
        function = inverse(obj, X)
            
        end

        function = constraint(obj, X)

        end
        function = act(obj, x)

        end

        function tau = vee(obj, tau_hat)

        end
        function tau_hat = hat(obj, tau)

        end

        function X = exp(obj, tau_hat)

        end
        function tau_hat = log(obj, X)

        end

        function = obj = cayley(obj, X)

        end
        function obj = inv_cayley(obj, X)

        end
    end
end

