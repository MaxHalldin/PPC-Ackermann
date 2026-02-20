function [dydt, inputs, APC] = ydotPPC2(t, Y, ref, params, rhoparams)
%   ydot outputs the first time derivative of the flat Ackermann system state
%   vector. Also outputs the motor force and rate of change of steering angle.
%   t       - 1x1 time. Used for trajectory functions.
%   Y       - 6x1 State vector of flat system. 
%   ref     - 6x1 cell array containing the reference trajectory functions.
%   params  - 7x1 cell array of vehicle parameters. (length & mass of axles)
%   rhos    - 3x1 cell array containing performance functions.
arguments (Input)
    t
    Y
    ref
    params
    rhoparams
end

arguments (Output)
    dydt
    inputs
    APC
end
    [L, M1, M2, k, lambda, ubar, lu, lo] = params{:};
    [rho0, rhoinf] = rhoparams{:};
    R = cellfun(@(c) c(t), ref);
    r = @(z) (rho0 - rhoinf)*exp(-lo*z) + rhoinf;

    A = [zeros(4,2) eye(4,4); zeros(2,6)];
    B = [zeros(4,2); eye(2,2)];

    %% Following section should probabaly be moved to a mother function or 
    % just outside this funciton
    e = Y(1) - R(1);

    a = poly(lambda);
    assert(r(0) > abs(a*e));
    rhoinf = einf*prod(lambda);
    assert(lu < min(lambda))

    %%
    Gamma = gammmas./(r.^(1:3)).';
    assert(k>0)
    xi      = a*sat(e,ebar)/r(t);
    ud      = -k*Jf(xi)*Tf(xi);
    u       = sat(ud, ubar);
    rhodot  = -lu*(r/(t) - rhoinf) + (u-ud)/xi;
    ehatdot = A*ehat + Gamma*Tf(eta(e,ehat,r(t)));

    APC = [rhodot ehatdot]; % Also needs to be integrated to create the 

    dydt = A*Y+B*u;

    %%
    u1 = sqrt(Y(3)^2+Y(4)^2);
    phi = atan((Y(3)*Y(6)-Y(4)*Y(5))/(u1^3)*L);
    u1dot = (Y(3)*Y(5)+Y(4)*Y(6))/u1;

    phidot = (L*sqrt(Y(3)^2 + Y(4)^2) * (u(2) * Y(3) - Y(4) * u(1)) *...
             (Y(3)^2 + Y(4)^2) - 3*(Y(5) * Y(3) + Y(6) * Y(4)) *...
             (Y(6) * Y(3) - Y(5) * Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + L^2 *...
             (Y(6) * Y(3) - Y(5)*Y(4))^2);
    F = (M1+M2+M2*tan(phi)^2)*u1dot + M2*phidot*sin(phi)*u1/(cos(phi)^3);

    inputs = [phidot; F];
    % disp("hit")
end

function res = Tf(x) % x in (-1,1)
    res = zeros(length(x),1);
    for i=1:length(x)
        assert((x(i) > -1 && x(i) < 1), "Argument %d is out of bounds", i);
        res(i) = 0.5*log( (1+x(i) / (1-x(i))) );
    end
end

function res = Jf(x) % x in (1,inf)
    res = zeros(length(x),1);
    for i=1:length(x)
        assert((x(i) > -1 && x(i) < 1), "Argument %d is out of bounds", i);
        res(i) = 1/(1-x(i)^2);
    end
end

function res = eta(e,ehat,r)
    res = zeros(length(e),1);
    for i=1:length(x)
        res(i) = (e(i)-ehat)/r;
    end
end

function res = sat(x, xbar)
    res = zeros(length(x),1);
    for i=1:length(x)
        res(i) = sign(x(i))*min([x(i) xbar]);
    end
end