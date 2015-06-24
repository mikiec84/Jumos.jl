# Copyright (c) Guillaume Fraux 2014-2015
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# ============================================================================ #
#             Universe type, holding all the simulated system data
# ============================================================================ #
export Universe
export set_cell!, add_interaction!, get_masses!, check_masses
export set_positions!, set_velocities!

include("Topology.jl")
include("UnitCell.jl")
include("Interactions.jl")

type Universe
    cell::UnitCell
    topology::Topology
    interactions::Interactions
    positions::Array3D{Float64}
    velocities::Array3D{Float64}
    masses::Vector{Float64}
    data::Dict{Symbol, Any}
end

function Universe(cell::UnitCell, topology::Topology)
    universe = Universe(cell,
                    topology,
                    Interactions(),
                    Array3D(Float64, 0),
                    Array3D(Float64, 0),
                    Float64[],
                    Dict{Symbol, Any}()
            )
    universe.data[:universeerse] = universe
    universe.data[:step] = 0
    return universe
end

function Base.size(u::Universe)
    assert(size(u.positions, 2) == size(u.topology))
    return size(u.positions, 2)
end

@doc "
`get_masses!(universeerse)`: get masses from the topology, and store them in the
universeerse internal data. Returns the masses array.
" ->
function get_masses!(u::Universe)
    u.masses = atomic_masses(u.topology)
    return u.masses
end

@doc "
`check_masses(universeerse)`: Check that all masses are defined and are not equals to 0.
" ->
function check_masses(universe::Universe)
    if countnz(universe.masses) != size(universe.topology)
        bad_masses = Set()
        for (i, val) in enumerate(universe.masses)
            if val == 0.0
                union!(bad_masses, [universe.topology[i].name])
            end
        end
        missing = join(bad_masses, " ")
        throw(JumosError(
                "Missing masses for the following atomic types: $missing."
            ))
    end
end

function add_atom!(u::Universe, atom::Atom)
    add_atom!(u.topology, atom)
end

function add_liaison!(u::Universe, atom_i::Atom, atom_j::Atom)
    add_liaison!(u.topology, atom_i, atom_j)
end

function remove_atom!(u::Universe, index)
    remove_atom!(u.topology, index)
end

function remove_liaison!(u::Universe, atom_i::Atom, atom_j::Atom)
    remove_liaison!(u.topology, atom_i, atom_j)
end

@doc "
`set_positions!(universeerse)`: get masses from the topology, and store them in the
universeerse internal data. Returns the masses array.
" ->
function set_positions!(universeerse::Universe, positions::Array3D{Float64})
    universeerse.positions = positions
end

@doc "
`set_cell!(universeerse, cell)`

Set the universeerse unit cell. `cell` can be an UnitCell instance or a list of cell
lenghts and angles.

`set_cell!(universeerse, celltype, params)`

Set the universeerse unit cell to a cell with type `celltype` and cell parameters from
`params`
" ->
function set_cell!(universe::Universe, cell::UnitCell)
    sim.cell = cell
end

function set_cell!(universe::Universe, params)
    return setcell!(sim, UnitCell(params...))
end

function set_cell!{T<:Type{AbstractCellType}}(universe::Universe, celltype::T, params = tuple())
    return setcell!(sim, UnitCell(celltype, params...))
end

# Todo: Way to add a catchall interaction
function add_interaction!(universe::Universe, pot::PotentialFunction, atoms...;
                         computation=:auto, kwargs...)
    atoms_id = get_atoms_id(universe.topology, atoms...)
    computation = get_computation(pot; computation=computation, kwargs...)
    push!(universe.interactions, computation, atoms_id...)
    return nothing
end

include("Distances.jl")
include("chemharp.jl")
