function [res,U,inputs,s] = URED_APC(t,y,ref,params)
%URED_APC One step of URED APC for Ackermann vehicle
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
    
    % URED
    sigmax = ehat(1)-e(1);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(1:2:4)-params.k.*phix).';
    
    sigmay = ehat(2)-e(2);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(2:2:4)-params.k.*phiy).';
    ehatdot1 = [ehatdotx;ehatdoty];
    ehatdot1 = ehatdot1(:);
    
    % extended observer.
    sigmax = -params.dt*ehatdot1(3);
    phix = [ abs(sigmax)^(1/2)*sign(sigmax) + params.mu*abs(sigmax)^(3/2)*sign(sigmax);
            0.5*sign(sigmax) + 2*params.mu*sigmax + 1.5*params.mu^2*abs(sigmax)^2*sign(sigmax) ];
    ehatdotx = ([0 1; 0 0]*ehat(3:2:6)-params.k.*phix).';
    
    sigmay = -params.dt*ehatdot1(4);
    phiy = [ abs(sigmay)^(1/2)*sign(sigmay) + params.mu*abs(sigmay)^(3/2)*sign(sigmay);
            0.5*sign(sigmay) + 2*params.mu*sigmay + 1.5*params.mu^2*abs(sigmay)^2*sign(sigmay) ];
    ehatdoty = ([0 1; 0 0]*ehat(4:2:6)-params.k.*phiy).';
    ehatdot2 = [ehatdotx;ehatdoty];
    ehatdot2 = ehatdot2(:);
    
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
    u1  = sqrt(Y(3)^2+Y(4)^2)+eps;
    wd  = (Y(3)*Y(5)+Y(4)*Y(6))/(u1+eps);
    phi = atan(params.l * (Y(3)*Y(6)-Y(4)*Y(5)) / u1^3+eps );
    
    phidotd = (params.l*sqrt(Y(3)^2 + Y(4)^2) * (ud(2)*Y(3) - Y(4)*ud(1))...
        * (Y(3)^2 + Y(4)^2) - 3*(Y(5)*Y(3) + Y(6)*Y(4)) * (Y(6)*Y(3)...
        - Y(5)*Y(4))) / ((Y(3)^2 + Y(4)^2)^3 + params.l^2 * (Y(6) * Y(3)...
        - Y(5)*Y(4))^2+eps);
    
    % Unconstrained jerk
    wdotd = 1/(eps+sqrt(Y(3)^2+Y(4)^2))^3 * ...
        ( (Y(3)^2+Y(4)^2) * (Y(5)^2+Y(6)^2 + ud(1)*Y(4)+ud(2)*Y(3)) - ...
        ( Y(3)*Y(5) + Y(4)*Y(6) )^2 );
    
    
    Fd = (params.M1+params.M2+params.M2*tan(phi)^2)*wd + ...
        params.M2*phidotd*sin(phi)*u1/(cos(phi)^3+eps);
    
    % Input Saturation
    phidot = sat(phidotd, params.phidotbar);
    F = sat(Fd, params.Fbar);
    w = (F-u1*params.M2*phidot*sin(phi)/(cos(phi)^3)) / ...
        (params.M1 + params.M2 + params.M1*tan(phi));
    
    a = u1*wdotd+u1^2*w^2-w^2-u1^4*tan(phi)^2/params.l^2;
    b = (u1/params.l)*( phidotd*u1*(1+tan(phi)^2) + u1^2*w*tan(phi) );
    
    u(1) = (a-b)/(2*Y(4)+eps);
    u(2) = (a+b)/(2*Y(3)+eps);
    u = u.';
    
    
    cond = ( (abs(u(1))<=abs(ud(1))) && (abs(u(2))<=abs(ud(2))) );
    u = cond*u + ~cond*ud; % Use ud if inputs are not saturated.
    U = [u; ud];
    
    inputs = [F; Fd; phidot; phidotd];
    
    rhodot = (u-ud)./xi - params.lu.*(rho - params.rho_inf);
    
    ydot = params.A*Y + params.B*u;
    
    res = [ydot; ehatdot; rhodot];
end
