function s = INTEGRATION(s, dsolution_dt, w)
% INTEGRATION  Apply a weighted Euler update to the conserved variables.
%
%   s = INTEGRATION(s, dsolution_dt, w)
%
%   Advances each conserved field (rho, rho_u, rho_v, rho_E, and
%   optionally gamma_star and cv_star for non-equilibrium chemistry) by
%   the increment  w * flux.  When shock fitting is active the fluxes
%   are masked by s.shock.flow_cells and the shock position is updated
%   via MOVE_SHOCK.
%
%   Inputs:
%       s     - struct : Current flow s with conservative
%                                variables and grid / shock data.
%       dsolution_dt - struct : Time-derivative (flux) data containing
%                                .flux.rho, .flux.rho_u, .flux.rho_v,
%                                .flux.rho_E, and (when applicable)
%                                .flux.gamma_star, .flux.cv_star.
%       w            - double : Integration weight (e.g. dt for forward
%                                Euler, or a Runge-Kutta stage weight).
%
%   Outputs:
%       s     - struct : Updated s with incremented
%                                conservative variables (interior cells
%                                only, indices 2:end-1).
%
%   Notes:
%       - For shock-fitted runs (s.shock.enabled == true), fluxes are
%         multiplied by s.shock.flow_cells to zero out non-flow regions,
%         and MOVE_SHOCK is called to update the shock geometry.
%       - Non-equilibrium chemistry fields (gamma_star, cv_star) are
%         integrated only when s.chemistry.chemical_equilibrium is false.
%
% Part of: Hypersonics Stability MATLAB Solver - Time Marching Module

    if s.shock.enabled
        %% Shock-fitted integration (masked by flow_cells)
        s.var.rho(2:end-1,2:end-1)   = s.var.rho(2:end-1,2:end-1)   + w * dsolution_dt.flux.rho   .* s.shock.flow_cells;
        s.var.rho_u(2:end-1,2:end-1) = s.var.rho_u(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_u .* s.shock.flow_cells;
        s.var.rho_v(2:end-1,2:end-1) = s.var.rho_v(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_v .* s.shock.flow_cells;
        s.var.rho_E(2:end-1,2:end-1) = s.var.rho_E(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_E .* s.shock.flow_cells;

        if ~s.chemistry.chemical_equilibrium
            s.var.gamma_star(2:end-1,2:end-1) = s.var.gamma_star(2:end-1,2:end-1) + w * dsolution_dt.flux.gamma_star .* s.shock.flow_cells;
            s.var.cv_star(2:end-1,2:end-1)    = s.var.cv_star(2:end-1,2:end-1)    + w * dsolution_dt.flux.cv_star    .* s.shock.flow_cells;
        end

        % Update shock geometry
        s = MOVE_SHOCK(s, dsolution_dt, w);
    else
        %% Standard integration (no shock fitting)
        s.var.rho(2:end-1,2:end-1)   = s.var.rho(2:end-1,2:end-1)   + w * dsolution_dt.flux.rho;
        s.var.rho_u(2:end-1,2:end-1) = s.var.rho_u(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_u;
        s.var.rho_v(2:end-1,2:end-1) = s.var.rho_v(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_v;
        s.var.rho_E(2:end-1,2:end-1) = s.var.rho_E(2:end-1,2:end-1) + w * dsolution_dt.flux.rho_E;

        if ~s.chemistry.chemical_equilibrium
            s.var.gamma_star(2:end-1,2:end-1) = s.var.gamma_star(2:end-1,2:end-1) + w * dsolution_dt.flux.gamma_star;
            s.var.cv_star(2:end-1,2:end-1)    = s.var.cv_star(2:end-1,2:end-1)    + w * dsolution_dt.flux.cv_star;
        end
    end

end