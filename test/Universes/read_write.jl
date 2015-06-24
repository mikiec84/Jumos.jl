facts("Read and write Universes to files") do
    tmpname = tempname() * ".xyz"
    content = """6

    He 0 0 0
    He 10 30 56
    He 0 0 0
    H 0 0 0
    O 0 0 0
    Zn 0 0 0
    """
    open(tmpname, "w") do file
        write(file, content)
    end

    context("Read") do
        universe = Universe(tmpname)
        @fact size(universe) => 6
        @fact size(universe.topology) => 6

        @fact universe.topology[1].label => :He
        @fact universe.topology[6].label => :Zn

        @fact universe.positions[1] => [0, 0, 0]
        @fact universe.positions[2] => [10, 30, 56]
    end
    rm(tmpname)
end
