"""
    GET_T_REF(s)

Compute a reference advection time scale for energy normalization.

Estimates a characteristic advection time by dividing the total
post-shock mass by the mass flux through the shock. This provides a
non-temporal scaling factor for the transient growth gain when a
time-independent normalization is preferred.

# Arguments
- `s`: Solution Dict{String,Any} containing:
  - `s["mesh"]["volume"]` - Cell volumes
  - `s["var"]["rho"]` - Density field (with ghost cells)
  - `s["shock"]["flow_cells"]` - Active-cell mask
  - `s["mesh"]["Nchi"]` - Number of streamwise cells
  - `s["mesh"]["bt_area"]` - Cell boundary areas
  - `s["mesh"]["bt_y_normal"]` - Wall-normal component of boundary normals
  - `s["shock"]["cell_indices"]` - Indices of shocked cells

# Returns
- `T_ref`: Reference advection time (mass_downstream / inflow_mass_flux).

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function GET_T_REF(s::Dict{String,Any})
    ## Compute total post-shock mass
    volume = s["mesh"]["volume"]
    rho = s["var"]["rho"][2:end-1, 2:end-1]
    flow_cells = s["shock"]["flow_cells"]
    mass_downstream = sum(rho .* flow_cells .* volume)

    ## Compute inflow mass flux through the shock
    inflow_mass_flux = 0.0
    for i in 1:s["mesh"]["Nchi"]
        inflow_mass_flux += s["mesh"]["bt_area"][i, s["shock"]["cell_indices"][i, 1]] * s["mesh"]["bt_y_normal"][i, s["shock"]["cell_indices"][i, 1]]
    end

    ## Reference time = mass / flux
    T_ref = mass_downstream / inflow_mass_flux
    return T_ref
end
