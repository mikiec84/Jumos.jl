# Copyright (c) Guillaume Fraux 2014-2015
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# ============================================================================ #
#                       Compute interesting values
# ============================================================================ #

using Jumos.Constants #kB

# abstract Compute -> Defined in Simulations.jl
export Compute
export TemperatureCompute, VolumeCompute, EnergyCompute
#export PressureCompute

function have_compute{T<:Compute}(sim::Simulation, compute_type::Type{T})
	for compute in sim.computes
		if isa(compute, compute_type)
			return true
		end
	end
	return false
end

# ============================================================================ #

@doc "
Compute the temperature of a simulation using the relation
	T = 1/kB * 2K/(3N) with K the kinetic energy
" ->
immutable TemperatureCompute <: Compute end

function Base.call(::TemperatureCompute, universe::Universe)
	T = 0.0
    K = kinetic_energy(universe)
    natoms = size(universe)
	T = 1/kB * 2*K/(3*natoms)
	universe.data[:temperature] = T
	return T
end

# ============================================================================ #

#@doc "
#Compute the pressure of the system.
#" ->
#immutable PressureCompute <: Compute end

# ============================================================================ #

@doc "
Compute the volume of the current simulation cell
" ->
immutable VolumeCompute <: Compute end

function Base.call(::VolumeCompute, universe::Universe)
    V = volume(universe.cell)
	universe.data[:volume] = V
    return V
end

# ============================================================================ #

@doc "
Compute the energy of a simulation.
    EnergyCompute()(simulation::MolecularDynamic) returns a tuple
    (Kinetic_energy, Potential_energy, Total_energy)
" ->
immutable EnergyCompute <: Compute end

function Base.call(::EnergyCompute, universe::Universe)
    K = kinetic_energy(universe)
    P = potential_energy(universe)
	universe.data[:E_kinetic] = K
	universe.data[:E_potential] = P
	universe.data[:E_total] = P + K
	return K, P, P + K
end

function kinetic_energy(universe::Universe)
    K = 0.0
	const natoms = size(universe)
	@inbounds for i=1:natoms
		K += 0.5 * universe.masses[i] * norm2(universe.velocities[i])
	end
    return K
end

function potential_energy(universe::Universe)
    E = 0.0
	const natoms = size(universe)
	@inbounds for i=1:natoms, j=(i+1):natoms
        atom_i = universe.topology.atoms[i]
        atom_j = universe.topology.atoms[j]
		for potential in pairs(universe.interactions, atom_i, atom_j)
			E += potential(distance(universe, i, j))
		end
		# TODO: bonds, angles, ...
	end
    return E
end
