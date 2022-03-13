function x = numerical_fitting(point)

%% Approximation parameters for the auxiliary functions
alpha_phi = 1e-4;
alpha_l = 1e-12;
alpha_psi = 0.1;
beta_psi = 0.5;
alpha_h = 1e-7;

%% auxiliary functions

    function value = phi(x, y, t, alpha)
        u = max(x,y);
        v = min(x,y);
        
        if t > alpha
            value = u .* g(x,y,t);
        else
            value = u .* ((g(x,y,alpha)-1).*t/alpha + 1);
        end
        assert(sum(isnan(value)) == 0) 
    end

    function value = g(x, y, t)
        idx = (x == 0) & (y == 0);
        u = max(x,y);
        v = min(x,y);
        value = (1 + (v./u).^(1/t)).^t;
        value(1,idx) = 0;
        assert(sum(isnan(value)) == 0)
    end

    function value = l(t, alpha)
        value = zeros(size(t));
        idx = find(t >= alpha);
        value(1, idx) = t(1, idx) .* log(t(1, idx));
        assert(sum(isnan(value)) == 0)
    end

    function value = f(r, t, alpha, beta)
        value = 1 ./ (1 + r.^(1/t));
        if t < alpha
            idx = find((1-r) < beta);
            value(1,idx) = 1 ./ (1 + r(1,idx).^(1/alpha));
        end
        assert(sum(isnan(value)) == 0)
    end        
 
    function value = psi(x, y, t, alpha, beta)
        value = nan(size(x));
        value(1,(x == 0) & (y == 0)) = 1/2;
        
        idx1 = find((x-y) >= 0 & x ~=0);
        x1 = x(1, idx1);
        y1 = y(1, idx1);
        rr = y1./x1;
        value(1, idx1) = f(rr, t, alpha, beta);
        
        idx2 = find((y-x) > 0 & y ~= 0);
        x2 = x(1, idx2);
        y2 = y(1, idx2);
        rr = x2./y2;
        value(1, idx2) = 1 - f(rr, t, alpha, beta);
        assert(sum(isnan(value)) == 0)
    end

    function value = H(x,y,z, a1, a2, a3, e1, e2, alpha)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = phi(phi(x_bar, y_bar, e2, alpha), z_bar, e1, alpha);
        assert(sum(isnan(value)) == 0)
    end

    function value = F(x,y,z, a1, a2, a3, e1, e2, alpha, HH)
        value = sqrt(a1*a2*a3) * (HH -1);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = -2/a1 * HH .*  ...
            psi(phi(x_bar, y_bar, e2, alpha_phi), z_bar, e1, alpha_psi, beta_psi) .* psi(x_bar, y_bar, e2, alpha_psi, beta_psi);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = -2/a2 * HH .*  ...
            psi(phi(x_bar, y_bar, e2, alpha_phi), z_bar, e1, alpha_psi, beta_psi) .* psi(y_bar, x_bar, e2, alpha_psi, beta_psi);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = -2/a3 * HH .* ...
            psi(z_bar, phi(x_bar, y_bar, e2, alpha_phi), e1, alpha_psi, beta_psi); 
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_de1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = -(HH .* ...
            (l(psi(z_bar, phi(x_bar, y_bar, e2, alpha_phi), e1, alpha_psi, beta_psi), alpha_l) + ...
             l(psi(phi(x_bar, y_bar, e2, alpha_phi), z_bar, e1, alpha_psi, beta_psi), alpha_l)));
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_de2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH)
        x_bar = (x/a1).^2;
        y_bar = (y/a2).^2;
        z_bar = (z/a3).^2;
        
        value = (-HH .* ...
            psi(phi(x_bar, y_bar, e2, alpha_phi), z_bar, e1, alpha_psi, beta_psi) .* ...
            (l(psi(x_bar,y_bar,e2, alpha_psi, beta_psi), alpha_l) + l(psi(y_bar,x_bar,e2, alpha_psi, beta_psi), alpha_l)));
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dx(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH)
        value = -x / (alpha_h^2*a1) .* dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
        idx = abs(x/a1) > alpha_h;
        value(1,idx) = -a1./x(1,idx) .* dH_da1(x(1,idx),y(1,idx),z(1,idx),a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH(1,idx));
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dy(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH)
        value = -y / (alpha_h^2*a2) .* dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
        idx = abs(y/a2) > alpha_h;
        value(1,idx) = -a2./y(1,idx) .* dH_da2(x(1,idx),y(1,idx),z(1,idx),a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH(1,idx));
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dz(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH)
        value = -z / (alpha_h^2*a3) .* dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
        idx = abs(z/a3) > alpha_h;
        value(1,idx) = -a3./z(1,idx) .* dH_da3(x(1,idx),y(1,idx),z(1,idx),a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH(1,idx));
        assert(sum(isnan(value)) == 0)
    end

    function value = dT_dgamma(gamma,beta,alpha,tx,ty,tz,x,y,z)
        value(1,:) = cos(beta)*sin(gamma)*(tx - x) - cos(beta)*cos(gamma)*(ty - y);
        value(2,:) = (cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma))*(tx - x) + (cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta))*(ty - y);
        value(3,:) = - (cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma))*(tx - x) - (sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta))*(ty - y);
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dT_dbeta(gamma,beta,alpha,tx,ty,tz,x,y,z)
        value(1,:) = cos(beta)*(tz - z) + cos(gamma)*sin(beta)*(tx - x) + sin(beta)*sin(gamma)*(ty - y);
        value(2,:) = sin(alpha)*sin(beta)*(tz - z) - cos(beta)*cos(gamma)*sin(alpha)*(tx - x) - cos(beta)*sin(alpha)*sin(gamma)*(ty - y);
        value(3,:) = cos(alpha)*sin(beta)*(tz - z) - cos(alpha)*cos(beta)*cos(gamma)*(tx - x) - cos(alpha)*cos(beta)*sin(gamma)*(ty - y);
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dT_dalpha(gamma,beta,alpha,tx,ty,tz,x,y,z)
        value(1,:) = zeros(size(x));
        value(2,:) = (cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma))*(ty - y) - (sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta))*(tx - x) - cos(alpha)*cos(beta)*(tz - z);
        value(3,:) = (cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma))*(ty - y) - (cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta))*(tx - x) + cos(beta)*sin(alpha)*(tz - z);
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dT_dtx(gamma,beta,alpha,tx,ty,tz,x,y,z)
        value(1,:) = -cos(beta)*cos(gamma)*ones(size(x));
        value(2,:) = (cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta)) * ones(size(x));
        value(3,:) = (- sin(alpha)*sin(gamma) - cos(alpha)*cos(gamma)*sin(beta)) * ones(size(x));
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dT_dty(gamma,beta,alpha,tx,ty,tz,x,y,z) 
        value(1,:) = -cos(beta)*sin(gamma)*ones(size(x));
        value(2,:) = (- cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma)) * ones(size(x));
        value(3,:) = (cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma)) * ones(size(x));
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dT_dtz(gamma,beta,alpha,tx,ty,tz,x,y,z)
        value(1,:) = sin(beta) * ones(size(x));
        value(2,:) = -cos(beta)*sin(alpha) * ones(size(x));
        value(3,:) = -cos(alpha)*cos(beta) * ones(size(x));
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dH_dX(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH)
        value = [dH_dx(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH);
                 dH_dy(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH);
                 dH_dz(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH)];
        assert(sum(isnan(value),'all') == 0)
    end

    function value = dH_dgamma(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)
        DT = dT_dgamma(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dbeta(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)     
        DT = dT_dbeta(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dalpha(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)
        DT = dT_dalpha(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dtx(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)
        DT = dT_dtx(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dty(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)
        DT = dT_dty(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dH_dtz(gamma,beta,alpha,tx,ty,tz, x,y,z, DH)
        DT = dT_dtz(gamma,beta,alpha,tx,ty,tz,x,y,z);
        
        value = dot(DH, DT);
        assert(sum(isnan(value)) == 0)
    end

    function value = dG_de1(xc,yc,zc, a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH)
        value = sqrt(a1*a2*a3) * dH_de1(xc,yc,zc, a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH);
    end

    function value = dG_de2(xc,yc,zc, a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH)
        value = sqrt(a1*a2*a3) * dH_de2(xc,yc,zc, a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH);
    end

    function value = dG_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        value = 0.5*(a1*a2*a3)^(-0.5)*a2*a3 * (HH - 1) + ...
                (a1*a2*a3)^(0.5)*dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
    end

    function value = dG_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        value = 0.5*(a1*a2*a3)^(-0.5)*a1*a3 * (HH - 1) + ...
                (a1*a2*a3)^(0.5)*dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
    end

    function value = dG_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH)
        value = 0.5*(a1*a2*a3)^(-0.5)*a2*a3 * (HH - 1) + ...
                (a1*a2*a3)^(0.5)*dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH);
    end

    function value = dG_dgamma(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dgamma(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end

    function value = dG_dbeta(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dbeta(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end

    function value = dG_dalpha(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dalpha(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end

    
    function value = dG_dtx(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dtx(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end

    function value = dG_dty(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dty(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end
        
    function value = dG_dtz(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3)
        value = (a1*a2*a3)^(0.5) * dH_dtz(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH);
    end    

%% cost function

    function [C,J] = myFun(para, X, alpha_phi, alpha_psi, beta_psi, alpha_l, alpha_h)
        e1 = para(1);
        e2 = para(2);
        a1 = para(3);
        a2 = para(4);
        a3 = para(5);
        
        gamma = para(6);
        beta = para(7);
        alpha = para(8);
        tx = para(9);
        ty = para(10);
        tz = para(11);
        
        R = eul2rotm([gamma, beta, alpha]);
        t = [tx;ty;tz];
        
        xx = X(1,:);
        yy = X(2,:);
        zz = X(3,:);
        
        Xc = R \ (X-t);
        xc = Xc(1,:);
        yc = Xc(2,:);
        zc = Xc(3,:);
        
        HH = H(xc,yc,zc, a1, a2, a3, e1, e2, alpha_phi);
        
        x_bar = (xc/a1).^2;
        y_bar = (yc/a2).^2;
        z_bar = (zc/a3).^2;
        
        phi1 = phi(x_bar, y_bar, e2, alpha_phi);
        
        psi1 = psi(phi1, z_bar, e1, alpha_psi, beta_psi);
        psi2 = psi(x_bar, y_bar, e2, alpha_psi, beta_psi);
        psi3 = psi(y_bar, x_bar, e2, alpha_psi, beta_psi);
        psi4 = psi(z_bar, phi1, e1, alpha_psi, beta_psi);
        
        
        Ha1 = -2/a1 * HH .* psi1 .* psi2;
        Ha2 = -2/a2 * HH .* psi1 .* psi3;
        Ha3 = -2/a3 * HH .* psi4; 
        
        Hx =  -xc / (alpha_h^2*a1) .* Ha1;
        idx = abs(xc/a1) > alpha_h;
        Hx(1,idx) = -a1./xc(1,idx) .* Ha1(1,idx);
        assert(sum(isnan(Hx)) == 0)
        
        Hy = -yc / (alpha_h^2*a2) .* Ha2;
        idx = abs(yc/a2) > alpha_h;
        Hy(1,idx) = -a2./yc(1,idx) .* Ha2(1,idx);
        assert(sum(isnan(Hy)) == 0)
        
        Hz = -zc / (alpha_h^2*a3) .* Ha3;
        idx = abs(zc/a3) > alpha_h;
        Hz(1,idx) = -a3./zc(1,idx) .* Ha3(1,idx);
        assert(sum(isnan(Hz)) == 0)
        
        DH  = [Hx;Hy;Hz];
        assert(sum(isnan(DH),'all') == 0)
        
        l1 = l(psi4, alpha_l);
        l2 = l(psi1, alpha_l);
        l3 = l(psi2, alpha_l);
        l4 = l(psi3, alpha_l);
        
        He1 = -(HH .* (l1 +l2));
        assert(sum(isnan(He1)) == 0)
        
        He2 = (-HH .* psi1 .* (l3 + l4));
        assert(sum(isnan(He2)) == 0)
        
        
        
        
%         DH = dH_dX(xc,yc,zc, a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH);
        
        C = F(xc,yc,zc, a1, a2, a3, e1, e2, alpha_phi, HH)';
        
        J = [sqrt(a1*a2*a3) * He1;
             sqrt(a1*a2*a3) * He2;
             0.5*(a1*a2*a3)^(-0.5)*a2*a3 * (HH - 1) + (a1*a2*a3)^(0.5)*Ha1;
             0.5*(a1*a2*a3)^(-0.5)*a1*a3 * (HH - 1) + (a1*a2*a3)^(0.5)*Ha2;
             0.5*(a1*a2*a3)^(-0.5)*a2*a3 * (HH - 1) + (a1*a2*a3)^(0.5)*Ha3;
             dG_dgamma(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);
             dG_dbeta(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);
             dG_dalpha(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);
             dG_dtx(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);
             dG_dty(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);
             dG_dtz(gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3);];
         
         
        J = J';

    end

%% lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display','off', 'SpecifyObjectiveGradient', true,...
                       'MaxIterations', 50);

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
% w = ones(1, size(point,2));
% W = ones(3, size(point,2));

euler0 = [rotm2eul([eigenVector(:, 1), eigenVector(:, 2), cross(eigenVector(:, 1), eigenVector(:, 2))]);
          rotm2eul([eigenVector(:, 1), eigenVector(:, 3), cross(eigenVector(:, 1), eigenVector(:, 3))]);
          rotm2eul([eigenVector(:, 2), eigenVector(:, 3), cross(eigenVector(:, 2), eigenVector(:, 3))]);];
x0 = [ones(3, 2), ones(3, 1) * median(abs(point(1, :))), ones(3, 1) * median(abs(point(2, :))), ones(3, 1) * median(abs(point(3, :))), euler0, zeros(3, 3)];
% x0 = [ones(1, 2), median(abs(point(1, :))), median(abs(point(2, :))), median(abs(point(3, :)))];
upper = 4 * max(max(abs(point)));
% lb = [0.0 0.0 0.01 0.01 0.01];
% ub = [2.0 2.0 ones(1, 3) * upper];
lb = [0.0 0.0 0.01 0.01 0.01 -2*pi -2*pi -2*pi -ones(1, 3) * upper];
ub = [2.0 2.0 ones(1, 3) * upper  2*pi 2*pi 2*pi ones(1, 3) * upper];

cost_func = @(x) myFun(x, point, alpha_phi, alpha_psi, beta_psi, alpha_l, alpha_h);

x = zeros(1, 11);
residue = Inf * ones(1, 3);
for i = 1 : size(x0, 1)
    [x(i, :), residue(i)] = lsqnonlin(cost_func, x0(i, :),lb,ub, options);
end
[~, I] = min(residue);
x = x(I, :);

% x = lsqnonlin(cost_func, x0,lb,ub, options);
x(9 : 11) = x(9 : 11) + t0';





end