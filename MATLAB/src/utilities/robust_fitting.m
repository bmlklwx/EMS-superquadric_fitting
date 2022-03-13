function [x, w] = robust_fitting(point)
% arguname = argnames(g);
%% init
t0 = mean(point, 2);
point = point - t0;
CovM = point * point' ./ size(point, 2);
[EigenVector, EigenValue] = eig(CovM);
[eigenValue, I] = sort(diag(EigenValue));
eigenVector = EigenVector(:,I);

% if abs(eigenValue(1) - eigenValue(2)) < abs(eigenValue(2) - eigenValue(3))
%     euler0 = rotm2eul([eigenVector(:,1),eigenVector(:,2), cross(eigenVector(:,1),eigenVector(:,2))]);
% else
%     euler0 = rotm2eul([eigenVector(:,2),eigenVector(:,3), cross(eigenVector(:,2),eigenVector(:,3))]);
% end
% x0 = [1, 1, median(abs(point(1, :))), median(abs(point(2, :))), median(abs(point(3, :))), euler0, zeros(1, 3)];
w = ones(1, size(point,2));
W = ones(3, size(point,2));

euler0 = [rotm2eul([eigenVector(:, 1), eigenVector(:, 2), cross(eigenVector(:, 1), eigenVector(:, 2))]);
          rotm2eul([eigenVector(:, 1), eigenVector(:, 3), cross(eigenVector(:, 1), eigenVector(:, 3))]);
          rotm2eul([eigenVector(:, 2), eigenVector(:, 3), cross(eigenVector(:, 2), eigenVector(:, 3))]);];
      
X0 = [ones(3, 2), ones(3, 1) * median(abs(point(1, :))), ones(3, 1) * median(abs(point(2, :))), ones(3, 1) * median(abs(point(3, :))), euler0, zeros(3, 3)];
X = zeros(3,11);
for k = 1:3
    x0 = X0(k,:);
    
    %% hyper-parameter
    alpha = 0.7; % between 0.5 and 1 0.7
    h = floor(alpha*size(point,2)); % h > 11
    % h = 15;
    goodness = 1e-3; % stopping condition 1e-3
    maxIter = 30; %30

    %% optimoption
    upper = 2 * max(max(abs(point)));
    lb = [0.1 0.1 0.01 0.01 0.01 -2*pi -2*pi -2*pi -ones(1, 3) * upper];
    ub = [2.0 2.0 ones(1, 3) * upper  2*pi 2*pi 2*pi ones(1, 3) * upper];
    cost_func = @(x) dist(x, point, w);
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display','off', 'MaxIterations', 2);

    %% optimize
    [x] = lsqnonlin(cost_func, x0, lb, ub, options);

    %% robust
    i = 0;
    while i < maxIter
        i = i + 1;
        Ri = dist(x, point, ones(1,size(point,2)))';
        mean_residual = sum(Ri.^2 .* w) / sum(w);
        while 1
            w = 2 - 0.5*Ri.^2/mean_residual;
            w(Ri.^2 < 2*mean_residual) = 1;
            w(Ri.^2 > 4*mean_residual) = 0;

            T = sum(w>0);
            if T < alpha  * size(point,2)
                mean_residual = 2 * mean_residual;
            else
                break
            end
        end

        idx = randperm(size(point,2));
        idx = idx(1:h);
        sub_point = point(:,idx);
        sub_w = w(idx);

    %     d = 0;
    %     for j = 1:h
    %         d = d - sub_w(j)^2*(double(subs(G,arguname,[x,sub_point(:,j)']))\...
    %             double(subs(g,arguname,[x,sub_point(:,j)'])));
    %     end
    %     
    %     x_new = x + d';

        cost_func = @(x) dist(x, sub_point, sub_w);

        x_new = lsqnonlin(cost_func, x, lb, ub, options);

        residual_old = mean(dist(x, point, w).^2);
        residual_new = mean(dist(x_new, point, w).^2);

        if residual_new < residual_old
            x = x_new;
        end

        if min(residual_old,residual_new) < goodness
            break
        end

    
    end

    x(9:11) = x(9:11) + t0'; 
    X(k,:) = x;
    W(k,:) = w;
    cost(k) = min(residual_old,residual_new);
end

[~,i] = min(cost);
x = X(i,:);
w = W(i,:);
    
%%
    function [value] = dist(para, X, weight)
        % Solina and Bajcsy, 1990
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
        
        X_c = R' * X - R' * t';
        value = (para(3) * para(4) * para(5)) ^ (1/ 2) .* (((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ para(1) - 1);
        
        value = value .* weight;
        value = value';
        value(isnan(value)) = 0;
        value(isinf(value)) = 0;
    end
    

end