function [E, h, u_num] = solve_one_m_SBP7_upwind(m, xl, xr, r_star, t_end, CFL, doPlot)
%SOLVE_ONE_M_SBP7_UPWIND  2D acoustics with 7th-order upwind SBP + Projection BCs.
%
%   u = [p; v; w] on Ω = [xl,xr]×[xl,xr] (same interval for x and y)
%   Boundary conditions: v=0 on West/East, w=0 on South/North
%   Time integration: RK4 with dt = CFL*h
%
% Inputs:
%   m      grid points per dimension
%   xl,xr  domain endpoints
%   r_star (unused here; kept for compatibility)
%   t_end  final time
%   CFL    CFL number for dt = CFL*h
%   doPlot true/false
%
% Outputs:
%   E      placeholder (set to NaN unless you add an energy diagnostic)
%   h      grid spacing
%   u_num  final state vector [p; v; w]

    if nargin < 7 || isempty(doPlot), doPlot = true; end %#ok<NASGU>
    %#ok<INUSD> r_star

    % ------------------------------------------------------------
    % Grid
    % ------------------------------------------------------------
    h  = (xr - xl) / (m - 1);
    x  = linspace(xl, xr, m);
    y  = linspace(xl, xr, m);
    [X, Y] = meshgrid(x, y);
    Nxy = m * m;                 % number of grid points in 2D

    % ------------------------------------------------------------
    % 1D SBP norm H and inverse HI (7th-order upwind)
    % ------------------------------------------------------------
    H = speye(m);
    H(1:6,1:6) = [0.19087e5 / 0.60480e5 0 0 0 0 0; 
                  0 0.84199e5 / 0.60480e5 0 0 0 0; 
                  0 0 0.18869e5 / 0.30240e5 0 0 0; 
                  0 0 0 0.37621e5 / 0.30240e5 0 0; 
                  0 0 0 0 0.55031e5 / 0.60480e5 0; 
                  0 0 0 0 0 0.61343e5 / 0.60480e5];

    H(m-5:m, m-5:m) = fliplr(flipud(H(1:6,1:6)));
    H  = H * h;
    HI = spdiags(1 ./ diag(H), 0, m, m);   % faster than inv(H) since H is diagonal

    % ------------------------------------------------------------
    % 1D SBP Qp, Qm (upwind)
    % ------------------------------------------------------------
    Qp = (-1/105*diag(ones(m-3,1),-3) ...
          +1/10*diag(ones(m-2,1),-2) ...
          -3/5*diag(ones(m-1,1),-1) ...
          -1/4*diag(ones(m,1),0) ...
          +1*diag(ones(m-1,1),+1) ...
          -3/10*diag(ones(m-2,1),+2) ...
          +1/15*diag(ones(m-3,1),+3) ...
          -1/140*diag(ones(m-4,1),+4));

    Q_U = [-0.265e3 / 0.300272e6   0.1587945773e10 / 0.2432203200e10  -0.1926361e7 / 0.25737600e8   -0.84398989e8 / 0.810734400e9   0.48781961e8 / 0.4864406400e10  0.3429119e7 / 0.202683600e9;
           -0.1570125773e10 / 0.2432203200e10  -0.26517e5 /0.1501360e7  0.240029831e9 / 0.486440640e9   0.202934303e9 / 0.972881280e9   0.118207e6 / 0.13512240e8  -0.231357719e9 / 0.4864406400e10;
            0.1626361e7 / 0.25737600e8   -0.206937767e9 / 0.486440640e9  -0.61067e5 / 0.750680e6      0.49602727e8 / 0.81073440e8   -0.43783933e8 / 0.194576256e9   0.51815011e8 / 0.810734400e9;
            0.91418989e8 / 0.810734400e9  -0.53314099e8 / 0.194576256e9  -0.33094279e8 / 0.81073440e8  -0.18269e5 / 0.107240e6       0.440626231e9 / 0.486440640e9  -0.365711063e9 / 0.1621468800e10;
           -0.62551961e8 / 0.4864406400e10  0.799e3 / 0.35280e5          0.82588241e8 / 0.972881280e9  -0.279245719e9 / 0.486440640e9  -0.346583e6 / 0.1501360e7       0.2312302333e10 / 0.2432203200e10;
           -0.3375119e7 / 0.202683600e9    0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9  0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5];

    Qp(1:6,1:6) = Q_U;
    Qp(m-5:m, m-5:m) = flipud(fliplr(Q_U))';  % per your operator definition
    Qm = -Qp';

    e1 = sparse(m,1); e1(1) = 1;
    em = sparse(m,1); em(m) = 1;

    Dp = HI * (Qp - 0.5*(e1*e1') + 0.5*(em*em'));
    Dm = HI * (Qm - 0.5*(e1*e1') + 0.5*(em*em'));

    % ------------------------------------------------------------
    % 2D derivative operators via Kronecker products
    % ------------------------------------------------------------
    I  = speye(m);
    Dp_x = kron(Dp, I);
    Dm_x = kron(Dm, I);
    Dp_y = kron(I, Dp);
    Dm_y = kron(I, Dm);

    Z = sparse(Nxy, Nxy);

    % System matrices (block form for u=[p; v; w])
    Dx = [ Z,   Dp_x, Z;
           Dm_x, Z,   Z;
           Z,   Z,   Z ];

    Dy = [ Z,   Z,   Dp_y;
           Z,   Z,   Z;
           Dm_y, Z,   Z ];

    % ------------------------------------------------------------
    % Projection to enforce boundary conditions
    %   v=0 on W/E, w=0 on S/N
    % ------------------------------------------------------------
    H2   = kron(H, H);
    HI2  = kron(HI, HI);
    HIbar = kron(speye(3), HI2);         % 3Nxy x 3Nxy

    % Boundary extraction operators:
    Ew = kron(I, e1');   % West:  x=xl for all y  -> m x Nxy
    Ee = kron(I, em');   % East:  x=xr for all y
    Es = kron(e1', I);   % South: y=xl for all x
    En = kron(em', I);   % North: y=xr for all x

    e_v = [0 1 0];
    e_w = [0 0 1];

    L = [ kron(e_v, Ew);
          kron(e_v, Ee);
          kron(e_w, Es);
          kron(e_w, En) ];             % (4m) x (3Nxy)

    A = L * HIbar * L.';               % (4m) x (4m)
    Pproj = speye(3*Nxy) - HIbar * L.' * (A \ L);

    % Semi-discrete operator
    M = -Pproj * (Dx + Dy) * Pproj;

    % ------------------------------------------------------------
    % Initial condition
    % ------------------------------------------------------------
    p0 = exp(-100 * (X.^2 + Y.^2));
    v0 = zeros(m, m);
    w0 = zeros(m, m);

    u = [p0(:); v0(:); w0(:)];

    % ------------------------------------------------------------
    % Time integration (RK4)
    % ------------------------------------------------------------
    dt = CFL * h;
    if t_end <= 0 || dt <= 0
        nSteps = 0;
    else
        nSteps = max(1, round(t_end / dt));
        dt = t_end / nSteps;  % make final time hit exactly
    end

    f = @(q) M * q;

    for n = 1:nSteps
        k1 = f(u);
        k2 = f(u + 0.5*dt*k1);
        k3 = f(u + 0.5*dt*k2);
        k4 = f(u + dt*k3);
        u  = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    u_num = u;
    E = NaN;  % TODO: add an energy diagnostic if you want (see below)

    % ------------------------------------------------------------
    % Plotting (optional)
    % ------------------------------------------------------------
    if doPlot
        p = reshape(u(1:Nxy), m, m);

        figure;
        surf(X, Y, p);
        shading interp;
        title(sprintf('p(x,y,t), m=%d, t=%.3g', m, t_end));
        xlabel('x'); ylabel('y'); zlabel('p');

        figure;
        imagesc(x, y, p);
        axis xy; colorbar;
        title(sprintf('p(x,y,t), m=%d, t=%.3g', m, t_end));
        xlabel('x'); ylabel('y');
    end
end
