# Create a dummy universe with n atoms.
function testing_universe_from_size(n)
    side = round(n^(1/3)) # integer such as n <= side^3

    positions = Array3D(Float64, n)
    nplaced = 0
    for i=0:side-1, j=0:side-1, k=0:side-1
        nplaced += 1
        positions[nplaced] = [2.0*i, 2.0*j, 2.0*k]
        nplaced == n ? break : nothing
    end

    top = Topology()
    for i=1:n
        add_atom!(top, Atom("He"))
    end
    univ = Universe(UnitCell(side*2.0), top)
    set_positions!(univ, positions)

    add_interaction!(univ, LennardJones(unit_from(0.8, "kJ/mol"), unit_from(2.0, "A")), "He", "He")
    create_velocities!(univ, 300)
    return univ
end

include("potentials.jl")
include("UnitCell.jl")
