function [res,U,inputs,s] = RED_APC(t,y,ref,params)
%RED_APC One step of URED APC for Ackermann vehicle
%   t       - Time, unused
%   y       - Flat state vector 
%   ref     - Reference trajectory 
%   params  - Parameters
arguments (Input)
    t
    y
    ref
    params
end

arguments (Output)
    res
    U
    inputs
    s
end

    Y = y(1:6);
    ehat = y(7:12);
    rho = y(13:14);
   
    e = Y(1:2)-ref;
    de = e(1:2)-ehat(1:2);

    % RED
    vx = [1/params.k(3)*params.DGF(params.k(3)^2*de(1));
         2*params.DGF(params.k(3)^2*de(1))*params.DGFprime(params.k(3)^2*de(1))];
    ehatdotx = ([ehat(3); 0] + params.k(1:2).*vx).';
    
    vy = [1/params.k(3)*params.DGF(params.k(3)^2*de(2));
         2*params.DGF(params.k(3)^2*de(2))*params.DGFprime(params.k(3)^2*de(2))];
    ehatdoty = ([ehat(4); 0] + params.k(1:2).*vy).';
    
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);

    % Extending observer
    
    dedot = params.dt*ehatdot1(3:4);
    k = params.kprim;

    vx = [1/k(3)*params.DGF(k(3)^2*dedot(1));
         2*params.DGF(k(3)^2*dedot(1))*params.DGFprime(k(3)^2*dedot(1))];
    ehatdotx = ([ehat(5); 0] + k(1:2).*vx).';
    
    vy = [1/k(3)*params.DGF(k(3)^2*dedot(2));
         2*params.DGF(k(3)^2*dedot(2))*params.DGFprime(k(3)^2*dedot(2))];
    ehatdoty = ([ehat(6); 0] + k(1:2).*vy).';
    
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);

    % Adding it all together
    ehatdot = [ehatdot1(1:2); ehatdot2];

    % Feedback control
    % 1. Linear filter
    s(1) = params.a*sat([e(1); ehat([3 5])],params.ebar);
    s(2) = params.a*sat([e(2); ehat([4 6])],params.ebar);
    s=s.';

    % 2. APC
    xi = s./(rho+eps);
    ud = -params.K*Jf(xi).*Tf(xi);

    % Vehicle control signals.
    
    u1  = sqrt(Y(3)^2+Y(4)^2);
    wd  = (Y(3)*Y(5)+Y(4)*Y(6))/(u1+eps);
    phi = atan(params.l * (Y(3)*Y(6)-Y(4)*Y(5)) / (u1^3+eps) );

    phidotd = (params.l*sqrt(Y(3)^2 + Y(4)^2) * (ud(2)*Y(3) - Y(4)*ud(1))...
        * (Y(3)^2 + Y(4)^2) - 3*(Y(5)*Y(3) + Y(6)*Y(4)) * (Y(6)*Y(3)...
        - Y(5)*Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + params.l^2 * (Y(6) * Y(3)...
        - Y(5)*Y(4))^2);
    
    wdotd = 1/(sqrt(Y(3)^2+Y(4)^2))^3 * ...
        ( (Y(3)^2+Y(4)^2) * (Y(5)^2+Y(6)^2 + ud(1)*Y(4)+ud(2)*Y(3)) - ...
        ( Y(3)*Y(5) + Y(4)*Y(6) )^2 );


    Fd = (params.M1+params.M2+params.M2*tan(phi)^2)*wd + ...
        params.M2*phidotd*sin(phi)*u1/(cos(phi)^3);
    
    % Input Saturation
    phidot = sat(phidotd, params.phidotbar);
    F = sat(Fd, params.Fbar);
    w = (F-params.M2*phidot*sin(phi)/(cos(phi)^3)) / ...
        (params.M1 + params.M2 + params.M1*tan(phi));

    a = u1*wdotd+u1^2*w^2-w^2-u1^4*tan(phi)^2/params.l^2;
    b = (u1/params.l)*( phidotd*u1*(1+tan(phi)^2) + u1^2*w*tan(phi) );

    u(1) = (a-b)/(2*Y(4)+eps);
    u(2) = (a+b)/(2*Y(3)+eps);
    u = u.';
    
    cond = ( (abs(u(1))<abs(ud(1))) && (abs(u(2))<abs(ud(2))) );
    u = cond*u + ~cond*ud; % Use ud if inputs are not saturated.
    U = [u; ud];

    inputs = [F; Fd; phidot; phidotd];

    rhodot = (u-ud)./xi - params.lu.*(rho - params.rho_inf);

    ydot = params.A*Y + params.B*u;

    res = [ydot; ehatdot; rhodot];
end
