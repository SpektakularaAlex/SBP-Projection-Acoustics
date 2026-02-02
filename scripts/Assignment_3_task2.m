function [E, h, u_num] = solve_one_m_SBP7_upwind_2D(m, xl, xr, t_end, CFL, doPlot)
%SOLVE_ONE_M_SBP7_UPWIND_2D  2D acoustics with 7th-order upwind SBP + projection BCs.
%
% System (example form):
%   p_t + v_x + w_y + beta(x) p = 0
%   v_t + (1/rho) p_x           = 0
%   w_t + (1/rho) p_y           = 0
%
% BCs imposed by projection:
%   v = 0 on West/East,  w = 0 on South/North
%
% Inputs:
%   m      grid points per dimension
%   xl,xr  domain endpoints (same for x and y)
%   t_end  final time
%   CFL    dt = CFL*h (then adjusted to hit t_end exactly)
%   doPlot true/false
%
% Outputs:
%   E      placeholder energy (NaN by default)
%   h      grid spacing
%   u_num  final solution vector [p; v; w] (length 3*m^2)

    if nargin < 6 || isempty(doPlot), doPlot = true; end

    % ----------------------------
    % Grid
    % ----------------------------
    h  = (xr - xl) / (m - 1);
    x  = linspace(xl, xr, m);
    y  = linspace(xl, xr, m);
    [X, Y] = meshgrid(x, y);

    N = m*m;           % points per component
    I = speye(m);

    % ----------------------------
    % 1D 7th-order upwind SBP: H, HI, Qp/Qm, D+/D-
    % ----------------------------
    [H, HI, Dp, Dm] = sbp7_upwind_1d(m, h);

    % ----------------------------
    % 2D derivative operators via Kronecker products
    % vec(U) uses MATLAB column-major ordering.
    % ----------------------------
    Dp_x = kron(Dp, I);
    Dm_x = kron(Dm, I);
    Dp_y = kron(I, Dp);
    Dm_y = kron(I, Dm);

    ZN = sparse(N, N);

    Dx = [ ZN,  Dp_x, ZN;
           Dm_x, ZN,  ZN;
           ZN,  ZN,  ZN ];

    Dy = [ ZN,  ZN,  Dp_y;
           ZN,  ZN,  ZN;
           Dm_y, ZN,  ZN ];

    % ----------------------------
    % Heterogeneous coefficients (example)
    %   - rho changes for x>0
    %   - beta(x) as damping on p
    % ----------------------------
    right = (X > 0);

    rho = ones(m,m);
    rho(right) = 0.5;                % example material jump

    beta_x = double(x(:) > 0);       % 1D beta depending on x only
    beta_vec = kron(ones(m,1), beta_x);   % N×1 (x varies with column index)
    BetaN = spdiags(beta_vec, 0, N, N);
    D_beta = blkdiag(BetaN, ZN, ZN);

    % ----------------------------
    % Material matrix C_inv (block-diagonal)
    % ----------------------------
    rho_vec = rho(:);
    c_vec   = ones(N,1);  % c=1 everywhere here; replace if needed

    C11 = spdiags(rho_vec .* (c_vec.^2), 0, N, N);  % rho*c^2
    C22 = spdiags(1 ./ rho_vec,          0, N, N);  % 1/rho
    C_inv = blkdiag(C11, C22, C22);

    % Combined spatial operator inside projection
    D = Dx + Dy + D_beta;

    % ----------------------------
    % Matrix-free projection P (enforce BCs)
    % ----------------------------
    % HI is diagonal -> HI2 and HI_bar are diagonal too.
    dhi  = full(diag(HI));         % m×1
    dhi2 = kron(dhi, dhi);         % N×1
    dhi_bar = [dhi2; dhi2; dhi2];  % 3N×1 diagonal of HI_bar

    % Boundary extraction operators for column-major vec():
    % West/East: first/last column => (e1' ⊗ I) vec(U)
    % South/North: first/last row  => (I ⊗ e1') vec(U)
    e1 = sparse(m,1); e1(1) = 1;
    em = sparse(m,1); em(m) = 1;

    Ew = kron(e1.', I);    % m×N
    Ee = kron(em.', I);
    Es = kron(I, e1.');    % m×N
    En = kron(I, em.');

    e_v = [0 1 0];
    e_w = [0 0 1];

    L = [ kron(e_v, Ew);
          kron(e_v, Ee);
          kron(e_w, Es);
          kron(e_w, En) ];       % (4m) × (3N)

    % A = L * HI_bar * L' without forming HI_bar:
    LHI = L * spdiags(dhi_bar, 0, 3*N, 3*N);  % (4m)×(3N)
    A = full(LHI * L.');                      % (4m)×(4m), small

    % Apply projection: P(u) = u - HI_bar * L' * (A \ (L*u))
    function Pu = applyP(u)
        y  = A \ (L*u);           % (4m)×1
        z  = L.' * y;             % (3N)×1
        Pu = u - (dhi_bar .* z);  % apply diagonal HI_bar
    end

    % RHS: u_t = - P * C_inv * D * P * u
    function ut = rhs(u)
        Pu  = applyP(u);
        ut  = -applyP(C_inv * (D * Pu));
    end

    % ----------------------------
    % Initial condition
    % ----------------------------
    p0 = exp(-100 * (X.^2 + Y.^2));
    v0 = zeros(m,m);
    w0 = zeros(m,m);
    u  = [p0(:); v0(:); w0(:)];

    % ----------------------------
    % RK4 time stepping
    % ----------------------------
    dt = CFL * h;
    if t_end <= 0 || dt <= 0
        nSteps = 0;
    else
        nSteps = max(1, round(t_end / dt));
        dt = t_end / nSteps;   % hit t_end exactly
    end

    for n = 1:nSteps
        k1 = rhs(u);
        k2 = rhs(u + 0.5*dt*k1);
        k3 = rhs(u + 0.5*dt*k2);
        k4 = rhs(u + dt*k3);
        u  = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    u_num = u;
    E = NaN;

    % ----------------------------
    % Plotting
    % ----------------------------
    if doPlot
        p = reshape(u(1:N), m, m);

        figure;
        surf(X, Y, p);
        shading interp;
        xlabel('x'); ylabel('y'); zlabel('p');
        title(sprintf('p(x,y,t), m=%d, t=%.3g', m, t_end));

        figure;
        imagesc(x, y, p);
        axis xy equal tight; colorbar;
        xlabel('x'); ylabel('y');
        title(sprintf('p(x,y,t), m=%d, t=%.3g', m, t_end));
    end
end

% =====================================================================
% Helper: 1D SBP7 upwind operator builder
% =====================================================================
function [H, HI, Dp, Dm] = sbp7_upwind_1d(m, h)
    % H (diagonal with boundary weights)
    H = speye(m);
    H(1:6,1:6) = [0.19087e5 / 0.60480e5 0 0 0 0 0; 
                  0 0.84199e5 / 0.60480e5 0 0 0 0; 
                  0 0 0.18869e5 / 0.30240e5 0 0 0; 
                  0 0 0 0.37621e5 / 0.30240e5 0 0; 
                  0 0 0 0 0.55031e5 / 0.60480e5 0; 
                  0 0 0 0 0 0.61343e5 / 0.60480e5];
    H(m-5:m, m-5:m) = fliplr(flipud(H(1:6,1:6)));
    H  = H * h;

    HI = spdiags(1 ./ diag(H), 0, m, m);

    % Interior stencil Qp (sparse banded)
    e = ones(m,1);
    Qp = (-1/105)*spdiags(e, -3, m, m) ...
       + ( 1/10)*spdiags(e, -2, m, m) ...
       + (-3/5 )*spdiags(e, -1, m, m) ...
       + (-1/4 )*spdiags(e,  0, m, m) ...
       + ( 1   )*spdiags(e, +1, m, m) ...
       + (-3/10)*spdiags(e, +2, m, m) ...
       + ( 1/15)*spdiags(e, +3, m, m) ...
       + (-1/140)*spdiags(e, +4, m, m);

    % Boundary closure (given)
    Q_U = [-0.265e3 / 0.300272e6   0.1587945773e10 / 0.2432203200e10  -0.1926361e7 / 0.25737600e8   -0.84398989e8 / 0.810734400e9   0.48781961e8 / 0.4864406400e10  0.3429119e7 / 0.202683600e9;
           -0.1570125773e10 / 0.2432203200e10  -0.26517e5 /0.1501360e7  0.240029831e9 / 0.486440640e9   0.202934303e9 / 0.972881280e9   0.118207e6 / 0.13512240e8  -0.231357719e9 / 0.4864406400e10;
            0.1626361e7 / 0.25737600e8   -0.206937767e9 / 0.486440640e9  -0.61067e5 / 0.750680e6      0.49602727e8 / 0.81073440e8   -0.43783933e8 / 0.194576256e9   0.51815011e8 / 0.810734400e9;
            0.91418989e8 / 0.810734400e9  -0.53314099e8 / 0.194576256e9  -0.33094279e8 / 0.81073440e8  -0.18269e5 / 0.107240e6       0.440626231e9 / 0.486440640e9  -0.365711063e9 / 0.1621468800e10;
           -0.62551961e8 / 0.4864406400e10  0.799e3 / 0.35280e5          0.82588241e8 / 0.972881280e9  -0.279245719e9 / 0.486440640e9  -0.346583e6 / 0.1501360e7       0.2312302333e10 / 0.2432203200e10;
           -0.3375119e7 / 0.202683600e9    0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9  0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5];

    Qp = sparse(Qp);
    Qp(1:6,1:6) = Q_U;
    Qp(m-5:m, m-5:m) = flipud(fliplr(Q_U)).';

    Qm = -Qp.';

    e1 = sparse(m,1); e1(1) = 1;
    em = sparse(m,1); em(m) = 1;

    Dp = HI * (Qp - 0.5*(e1*e1.') + 0.5*(em*em.'));
    Dm = HI * (Qm - 0.5*(e1*e1.') + 0.5*(em*em.'));
end
