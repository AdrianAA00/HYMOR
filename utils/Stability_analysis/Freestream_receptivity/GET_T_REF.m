function T_ref = GET_T_REF(s)
% GET_T_REF  Compute a reference advection time scale for energy normalization.
%
%   T_ref = GET_T_REF(s)
%
%   Estimates a characteristic advection time by dividing the total
%   post-shock mass by the mass flux through the shock. This provides a
%   non-temporal scaling factor for the transient growth gain when a
%   time-independent normalization is preferred.
%
%   Inputs:
%       s - Solution structure containing:
%                    .mesh.volume         - Cell volumes
%                    .var.rho             - Density field (with ghost cells)
%                    .shock.flow_cells    - Active-cell mask
%                    .mesh.Nchi           - Number of streamwise cells
%                    .mesh.bt_area        - Cell boundary areas
%                    .mesh.bt_y_normal    - Wall-normal component of boundary normals
%                    .shock.cell_indices  - Indices of shocked cells
%
%   Outputs:
%       T_ref - Reference advection time (mass_downstream / inflow_mass_flux)
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Compute total post-shock mass
    volume = s.mesh.volume;
    rho = s.var.rho(2:end-1, 2:end-1);
    flow_cells = s.shock.flow_cells;
    mass_dowstream = sum(rho .* flow_cells .* volume, "all");

    %% Compute inflow mass flux through the shock
    inflow_mass_flux = 0;
    for i = 1:s.mesh.Nchi
        inflow_mass_flux = inflow_mass_flux + s.mesh.bt_area(i, s.shock.cell_indices(i, 1)) * s.mesh.bt_y_normal(i, s.shock.cell_indices(i, 1));
    end

    %% Reference time = mass / flux
    T_ref = mass_dowstream ./ inflow_mass_flux;
end
