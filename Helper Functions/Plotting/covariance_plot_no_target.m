function [] = covariance_plot_no_target(x_pred, data, P_pred, labels, title)
%COVARIANCE_PLOT Summary of this function goes here
%   Detailed explanation goes here

n = size(P_pred, 1);

tiledlayout(n, n, "TileSpacing","compact");

for i = 1:n % column
    for j = 1:n % row
        if i - 1 < j
            nexttile((j - 1) * n + i)

            if i == j
                histogram(data(i, :), ceil(sqrt(size(data, 2))), FaceColor=[120, 120, 120]/256, EdgeColor=[120, 120, 120]/256, Normalization="pdf"); hold on

                std_pred = sqrt(P_pred(i, i));
                x = 2 * 3 * std_pred * (-1:0.01:1);
                fx_theory = pdf('Normal', x, 0, std_pred); %theoretical normal probability density
                hold on; plot(x + x_pred(i),fx_theory,'r'); %plot computed theoretical PDF  

                std_data = std(data(i, :));
                x = 2 * 3 * std_data * (-1:0.01:1);
                fx_theory = pdf('Normal', x, 0, std_data); %theoretical normal probability density
                hold on; plot(x + mean(data(i, :)),fx_theory,'b', LineStyle="--"); %plot computed theoretical PDF  

                set(gca,'YColor','none')

                xlim(x_pred(i) + [-6 * std_data, 6 * std_data])

                if i == 1 && j == 1
                    legend("Samples", "Prediction Distribution 3\sigma Ellipse", "Sample Distribution 3\sigma Ellipse")
                end
            else
                scatter(squeeze(data(i, :)), squeeze(data(j, :)), 2, [120, 120, 120]/256, "filled"); hold on
                plot_cov_ellipse(P_pred([i, j], [i, j]), [x_pred(i), x_pred(j)], Color = "r"); hold on
                plot_cov_ellipse(cov(data([i,j], :)'), [mean(data(i, :)), mean(data(j, :))], Color = "b", LineStyle = "--"); hold off
            end
    
            if i ~= 1
                set(gca,'YColor','none')
            else
                ylabel(labels(j));
            end
            if j ~= n
                set(gca,'XColor','none')
            else
                xlabel(labels(i))
            end
            grid on
        end
    end
end

sgtitle(title)

end

function [] = plot_cov_ellipse(P, xbar, options)
    arguments
        P
        xbar
        options.Color = "r"
        options.LineStyle = "-"
    end

    [P_eigvec, P_eigval] = eigs(P);
    thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
    ellipse_3sigma_data = P_eigvec * [3 * sqrt(P_eigval(1, 1)) * cos(thetas); 3 * sqrt(P_eigval(2, 2)) * sin(thetas)];
    plot(xbar(1) + ellipse_3sigma_data(1, :), xbar(2) + ellipse_3sigma_data(2, :), Color = options.Color, LineStyle = options.LineStyle);
end