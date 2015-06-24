# Copyright (c) Guillaume Fraux 2014
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# ============================================================================ #
#                    Convertion from and to Chemharp types
# ============================================================================ #
import Chemharp

"Convert a Chemharp Atom to a Jumos Atom"
function Atom(other::Chemharp.Atom)
    atom = Atom(Chemharp.name(other))
    atom.mass = Chemharp.mass(other)
    atom.charge = Chemharp.charge(other)
    return atom
end

"Convert a Chemharp UnitCell to a Jumos UnitCell"
function UnitCell(other::Chemharp.UnitCell)
    a, b, c = Chemharp.lengths(other)
    α, β, γ = Chemharp.angles(other)
    chrp_type = Chemharp.cell_type(other)
    if chrp_type == Chemharp.INFINITE
        celltype = InfiniteCell
    elseif chrp_type == Chemharp.ORTHOROMBIC
        celltype = OrthorombicCell
    elseif chrp_type == Chemharp.TRICLINIC
        celltype = TriclinicCell
    end
    return UnitCell(a, b, c, α, β, γ, celltype)
end

"Convert a Chemharp Topology to a Jumos Topology"
function Topology(other::Chemharp.Topology)
    topology = Topology()
    for i=0:Chemharp.natoms(other)-1
        atom = Atom(Chemharp.Atom(other, i))
        add_atom!(topology, atom)
    end
    nbonds = Chemharp.nbonds(other)
    # TODO: get bonds
    return topology
end

"Read a file, and use the first frame to build an Universe."
function Universe(filename::AbstractString; guess_bonds=false)
    traj = Chemharp.Trajectory(filename)
    frame = read(traj)
    if guess_bonds
        frame.guess_topology(true)
    end
    return Universe(frame)
end

function Universe(frame::Chemharp.Frame)
    cell = UnitCell(Chemharp.UnitCell(frame))
    topology = Topology(Chemharp.Topology(frame))
    universe = Universe(cell, topology)

    data = Array(Float32, 3, size(frame))
    Chemharp.positions!(frame, data)
    universe.positions = Array3D{Float64}(data)
    universe.velocities = Array3D(Float64, 0)
    if Chemharp.has_velocities(frame)
        Chemharp.velocities!(frame, data)
        universe.velocities = Array3D{Float64}(data)
    end

    get_masses!(universe)
    return universe
end
