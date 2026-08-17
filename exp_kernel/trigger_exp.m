function g_val = trigger_exp(child, parent, parameters, time_scalar)
% TRIGGER_EXP  Exponential triggering kernel.
%
%   g = K0 * strain * omega * exp(-omega * dt) * (rho/pi) * exp(-rho * r^2)
%
% Parameters: [lambda0, K0, omega, rho]
%   lambda0 : background rate (not used here)
%   K0      : excitation strength
%   omega   : temporal decay rate (higher = faster decay)
%   rho     : spatial decay rate (higher = faster decay)

K0    = parameters(2);
omega = parameters(3);
rho   = parameters(4);

% Preserve empty observation blocks when elapsed timestamps are available.
if isfield(child, 'time') && isfield(parent, 'time')
    dt = abs(child.time - parent.time);
else
    dt = abs(child.kk - parent.kk) * time_scalar;
end

% Spatial distance squared
r_sq = (child.firstx - parent.firstx)^2 + (child.firsty - parent.firsty)^2;

% Parent strain rate as mark
if parent.bb > 0
    strain = abs(parent.aa) / parent.bb;
else
    strain = 0;
    g_val = 0;
    return;
end

% Exponential kernel
temporal = omega * exp(-omega * dt);
spatial  = (rho / pi) * exp(-rho * r_sq);

g_val = K0 * strain * temporal * spatial;

end
