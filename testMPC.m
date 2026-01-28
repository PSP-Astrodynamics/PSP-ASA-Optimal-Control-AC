Isp_0 = 225;
I_0 =  [19150; 13600; 13600] * (1e-3) ^ 2;
Wind_0 = [0,1,0];


Isp = Isp_0 + randn(100,1).*(Isp_0*0.10);
I = I_0 + randn(100,1).*(I_0*0.10);
Wind = Wind_0 + randn(100,1).*(Wind_0*0.10);


for index = 0:length(Isp)
    run_mpc_func(Isp(index),I_0, [0,0,0], [0,0,0]);
end

for index = 0:length(Isp)
    run_mpc_func(Isp_0, I_0(index), [0,0,0], [0,0,0]);
end

for index = 0:length(Isp)
    run_mpc_func(Isp_0, I_0, Wind, [0,0,0]);
end