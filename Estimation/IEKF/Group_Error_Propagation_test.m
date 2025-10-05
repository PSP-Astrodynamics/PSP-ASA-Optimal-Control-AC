R = make_R2(deg2rad(90));
r = [0; -10];

X = SE2_RotationMatrix(R, r);

Sigma = diag([0.3, 0.01, 0.03]);

N = 10000;
r_array = zeros([2, N]);
for k = 1 : N
    X_sample(k) = sample_left_gaussian(Sigma, X);
    %X_sample(k) = lplus([1; 1e-5; 1e-5], X_sample(k));
    r_array(:, k) = X_sample(k).r;
end

P = Sigma(1:2, 1:2);
[P_eigvec, P_eigval] = eigs(P);
thetas = reshape(linspace(0, 2 * pi, N), 1, []);
ellipse_3sigma_data = P_eigvec * [3 * sqrt(P_eigval(1, 1)) * cos(thetas); 3 * sqrt(P_eigval(2, 2)) * sin(thetas)];
for k = 1 : N
    X_ellipse(k) = lplus([ellipse_3sigma_data(:, k); 1e-10], X);
    ellipse_3sigma_data(:, k) = X_ellipse(k).r;
end

[Xe,Ye,Ze] = ellipsoid(0,0,0,sqrt(Sigma(1,1)) * 3,sqrt(Sigma(2,2)) * 3,sqrt(Sigma(3,3)) * 3,40);

points = zeros([3, 41, 41]);
points(1, :, :) = Xe;
points(2, :, :) = Ye;
points(3, :, :) = Ze;

ellipsoid_3sigma_data = zeros([2, 41 * 41]);
for k = 1 : 41
    for j = 1 : 41
        X_ellipsoid = lplus([points(:, k, j) + 1e-5], X);
        ellipsoid_3sigma_data(:, (k - 1) * 41 + j) = X_ellipsoid.r;
    end
end

%%
scatter(r_array(1, :), r_array(2, :)); hold on
plot(ellipse_3sigma_data(1, :), ellipse_3sigma_data(2, :), LineWidth=3); hold on
plot(ellipsoid_3sigma_data(1, :), ellipsoid_3sigma_data(2, :), LineWidth=1.5)
grid on
axis equal

%%
R = make_R2(deg2rad(90));
r = [2; -10];

X = SE2_RotationMatrix(R, r);

Sigma = diag([0.1, 0.1, deg2rad(2)] .^ 2);

N = 10000;

r_array = zeros([2, N]);
for k = 1 : N
    X_sample(k) = sample_right_gaussian(X, Sigma);
    %X_sample(k) = lplus([1; 1e-5; 1e-5], X_sample(k));
    r_array(:, k) = X_sample(k).r;
end

P = Sigma(1:2, 1:2);
[P_eigvec, P_eigval] = eigs(P);
thetas = reshape(linspace(0, 2 * pi, N), 1, []);
ellipse_3sigma_data = P_eigvec * [3 * sqrt(P_eigval(1, 1)) * cos(thetas); 3 * sqrt(P_eigval(2, 2)) * sin(thetas)];
for k = 1 : N
    X_ellipse(k) = rplus(X, [ellipse_3sigma_data(:, k); 1e-10]);
    ellipse_3sigma_data(:, k) = X_ellipse(k).r;
end

[Xe,Ye,Ze] = ellipsoid(0,0,0,sqrt(Sigma(1,1)) * 3,sqrt(Sigma(2,2)) * 3,sqrt(Sigma(3,3)) * 3,40);

points = zeros([3, 41, 41]);
points(1, :, :) = Xe;
points(2, :, :) = Ye;
points(3, :, :) = Ze;

ellipsoid_3sigma_data = zeros([2, 41 * 41]);
for k = 1 : 41
    for j = 1 : 41
        X_ellipsoid = rplus(X, [points(:, k, j) + 1e-5]);
        ellipsoid_3sigma_data(:, (k - 1) * 41 + j) = X_ellipsoid.r;
    end
end

%%
scatter(r_array(1, :), r_array(2, :)); hold on
plot(ellipse_3sigma_data(1, :), ellipse_3sigma_data(2, :), LineWidth=3); hold on
plot(ellipsoid_3sigma_data(1, :), ellipsoid_3sigma_data(2, :), LineWidth=1.5)
grid on
axis equal

%%
for k = 1 : N
    tau = X_sample(k).Log;
    p(:, k) = X_sample(k).r;
    th(k) = rad2deg(tau(3));
end

histogram(th); hold on
xline(90 + rad2deg(sqrt(Sigma(3,3))))
xline(90 - rad2deg(sqrt(Sigma(3,3)))); hold off
