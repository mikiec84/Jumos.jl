# Copyright (c) Guillaume Fraux 2014-2015
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# ============================================================================ #
#                   Time integration takes place here
# ============================================================================ #

export Integrator

export VelocityVerlet, Verlet
# abstract Integrator -> Defined in MolecularDynamics.jl

function setup(::Integrator, ::Simulation, ::Universe) end

@doc "
Velocity Verlet integrator
" ->
type VelocityVerlet <: Integrator
    timestep::Float64
    accelerations::Array3D
end

function VelocityVerlet(timestep::Real)
    accelerations = Array3D(Float64, 0)
    return VelocityVerlet(timestep, accelerations)
end

function setup(integrator::VelocityVerlet, ::Simulation, universe::Universe)
    natoms = size(universe)
    if length(integrator.accelerations) != natoms
        resize!(integrator.accelerations, natoms)
        fill!(integrator.accelerations, 0)
    end
end

function Base.call(integrator::VelocityVerlet, md::MolecularDynamics, universe::Universe)
    const dt = integrator.timestep

    # Getting pointers to facilitate further reading
    positions = universe.positions
    velocities = universe.velocities
    accelerations = integrator.accelerations
    masses = universe.masses

    const natoms = size(universe)

    # Update velocities at t + ∆t/2
    @inbounds for i=1:natoms, dim=1:3
        velocities[dim, i] += 0.5*accelerations[dim, i]*dt
    end

    # Update positions at t + ∆t
    @inbounds for i=1:natoms, dim=1:3
        positions[dim, i] += velocities[dim, i]*dt
    end

    getforces!(md, universe)
    # Update accelerations at t + ∆t
    @inbounds for i=1:natoms, dim=1:3
        accelerations[dim, i] = md.forces[dim, i] / masses[i]
    end

    # Update velocities at t + ∆t
    @inbounds for i=1:natoms, dim=1:3
        velocities[dim, i] += 0.5*accelerations[dim, i]*dt
    end
end

# ============================================================================ #

@doc "
Basic Verlet integrator. Velocities are updated at t + 1/2 ∆t
" ->
type Verlet <: Integrator
    timestep::Float64
    tmp::Array3D     # Temporary array for computations
    prevpos::Array3D # Previous positions
    wrap_velocities::Bool
end

function Verlet(timestep::Real)
    return Verlet(timestep, Array3D(Float64, 0), Array3D(Float64, 0), false)
end

function setup(integrator::Verlet, sim::Simulation, universe::Universe)
    const natoms = size(universe)

    integrator.wrap_velocities = ispresent(sim, WrapParticles())

    if length(integrator.prevpos) != natoms || length(integrator.tmp) != natoms
        resize!(integrator.prevpos, natoms)
        resize!(integrator.tmp, natoms)

        dt = integrator.timestep
        # Approximate the positions at t - ∆t
        for i=1:natoms
            integrator.prevpos[i] = universe.positions[i] - universe.velocities[i].*dt
        end
    end
end

function Base.call(integrator::Verlet, md::MolecularDynamics, universe::Universe)
    const dt = integrator.timestep

    # Getting pointers to facilitate further reading
    positions = universe.positions
    velocities = universe.velocities
    prevpos = integrator.prevpos
    tmp = integrator.tmp
    masses = universe.masses

    natoms = size(universe)
    getforces!(md, universe)

    # Save positions at t
    @inbounds for i=1:natoms, dim=1:3
        tmp[dim, i] = positions[dim, i]
    end

    # Update positions at t + ∆t
    @inbounds for i=1:natoms, dim=1:3
        positions[dim, i] = (2.0 * positions[dim, i] - prevpos[dim, i] +
                                        (dt^2 / masses[i]) * md.forces[dim, i])
    end

    # Update velocities at t
    if integrator.wrap_velocities
        # If the postions are wrapped in the simulation, position is updated,
        # but not prevpos. So let's do it now.
        delta_pos = zeros(Float64, 3)
        @inbounds for i=1:natoms
            delta_pos = positions[i] - prevpos[i]
            minimal_image!(delta_pos, universe.cell)
            for dim=1:3
                velocities[dim, i] = delta_pos[dim] / (2.0 * dt)
            end
        end
    else
        @inbounds for i=1:natoms, dim=1:3
            velocities[dim, i] = (positions[dim, i] - prevpos[dim, i]) / (2.0 * dt)
        end
    end

    # Update saved position
    @inbounds for i=1:natoms, dim=1:3
        prevpos[dim, i] = tmp[dim, i]
    end
end
