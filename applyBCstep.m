function [p, rho, u, v, T] = applyBCstep(rho, u, v, T, R, T_inf, p, p_inf, u_inf, ...
                                      i_step, j_step, solidMask)
% applyBCs  Apply primitive‐variable boundary conditions, including flat bottom,
%           forward‐facing step, top freestream, left inflow, right outflow.
%
%  Inputs (same as before, plus):
%    i_step, j_step   – grid‐indices of the step “corner” (see text)
%    solidMask(nx,ny) – boolean mask: true = “solid” cell
%
%  Outputs: p, rho, u, v, T after BCs.

  [nx, ny] = size(rho);

  % ------------------------------------------------------
  %  3a) Bottom flat‐plate (only for i < i_step):
  %      Enforce no‐slip & isothermal at j = 1 for i=1..(i_step−1)
  for i = 1:(i_step-1)
    u(i,1) = 0;
    v(i,1) = 0;
    T(i,1) = T_inf;
    % Extrapolate pressure at wall (2nd order):
    if ny>=3
      p(i,1) = 2*p(i,2) - p(i,3);
    else
      p(i,1) = p_inf;
    end
    % Re‐compute rho from ideal‐gas law:
    rho(i,1) = p(i,1) / (R * T(i,1));
  end

  %  3b) Step top‐surface (for i = i_step..nx, at j = j_step):
  %      Also a no‐slip, isothermal wall
  for i = i_step:nx
    u(i,j_step) = 0;
    v(i,j_step) = 0;
    T(i,j_step) = T_inf;
    % Extrapolate pressure into the wall cell (2nd order from above):
    if (j_step+1 <= ny)
      p(i,j_step) = 2*p(i,j_step+1) - p(i,j_step+2);
    else
      p(i,j_step) = p_inf;
    end
    rho(i,j_step) = p(i,j_step)/(R * T(i,j_step));
  end

  % ----------------------------------------------
  %  3c) Mark all “solid” cells (below step) as exactly the same as the step wall.
  %      For i>=i_step, j<j_step, we simply do not update them.  It is convenient
  %      to re‐impose the wall BC at the cell centers immediately adjacent to that solid box.
  %
  %  In other words, if (i,j) is solidMask(i,j)==true, we skip updating it in the
  %  main MacCormack loops.  But to keep the primitive array “clean,” zero out velocities
  %  inside the solid, and set T = T_inf, rho/p to whatever you want (so they never blow up).
  for i = i_step:nx
    for j = 1:(j_step-1)
      if solidMask(i,j)
        u(i,j)   = 0;
        v(i,j)   = 0;
        T(i,j)   = T_inf;
        p(i,j)   = p_inf;          % or any large number, but p_inf is simplest
        rho(i,j) = p_inf./R./T_inf;        % maintain freestream density or p_inf/(R T_inf)
      end
    end
  end

  % ---------------------------------------------
  %  3d) Top freestream (j = ny):
  u(:,ny)   = u_inf;          
  v(:,ny)   = 0;
  T(:,ny)   = T_inf;
  p(:,ny)   = p_inf;
  rho(:,ny)= p(:,ny) ./ (R * T(:,ny));

  % ---------------------------------------------
  %  3e) Left inflow (i = 1):
  u(1,:)   = u_inf;          
  v(1,:)   = 0;
  T(1,:)   = T_inf;          
  p(1,:)   = p_inf;
  rho(1,:) = p(1,:) ./ (R * T(1,:));

  % ---------------------------------------------
  %  3f) Right outflow (i = nx), subsonic nonreflective:
  u(nx,:)   = 2*u(nx-1,:) - u(nx-2,:);
  v(nx,:)   = 2*v(nx-1,:) - v(nx-2,:);
  T(nx,:)   = 2*T(nx-1,:) - T(nx-2,:);
  p(nx,:)   = 2*p(nx-1,:) - p(nx-2,:);
  rho(nx,:) = p(nx,:) ./ (R .* T(nx,:));

end