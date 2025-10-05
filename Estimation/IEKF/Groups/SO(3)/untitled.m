g = sym("g", [3, 1]);

val = [zeros([3, 9]); 
       skew(g), zeros(3, 6); 
       zeros([3, 3]), eye(3), zeros([3, 3])]

expm(val)