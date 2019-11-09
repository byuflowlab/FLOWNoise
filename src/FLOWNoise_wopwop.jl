#=##############################################################################
# DESCRIPTION
    Tools for formation of PSU-WOPWOP's inputs/outputs.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################


"""
    `wopwop2vtk(grid, outputname, save_path; read_path="", paraview=false)`

Converts an output field from PSU-WOPWOP generated over a grid of observers to
vtk format. The user must make sure that the `grid` matches the exact grid that
WOPWOP used as an observer to generate the output fields.

# ARGUMENTS
* `grid::gt.AbstractGrid`        : Grid.
* `outputname::String`           : PSU-WOPWOP output file to read. For instance,
                                    `outputname=pressure` will look for the files
                                    pressure.nam and pressure.fn under path
                                    `read_path`.
* `save_path::String`           : Where to save vtk files.

# OPTIONAL ARGUMENTS
* `read_path::String`           : Path from where to read the inputs files.

"""
function wopwop2vtk(grid::gt.AbstractGrid, outputname::String, save_path::String;
                        read_path="", prompt=true,
                        verbose=true, v_lvl=0,
                        paraview=false)

    gt.create_path(save_path, prompt)

    if verbose
        println("\t"^v_lvl*"Reading output file...")
    end

    # Read header file
    f = open(joinpath(read_path, outputname*".nam"), "r")
    header = readlines(f)
    close(f)

    # Remove invalid characters from headers
    for str in [" ", "(", ")"]
        for i in 1:length(header)
            header[i] = replace(header[i], str, "")
        end
    end

    # Read function file
    f = open(joinpath(read_path, outputname*".fn"), "r")
    imax, jmax, tmax, fieldmax = [parse(elem) for elem in split(readline(f))]

        if verbose
        println("Found $(imax)x$(jmax) grid in WOPWOP output with"*
                    " $tmax time steps.")
    end

    # Read fields
    field = zeros(imax, jmax, tmax, fieldmax)
    for fieldi in 1:fieldmax # Iterate over fields

        if verbose
            println("\t"^(v_lvl+1)*"Reading field $(header[fieldi])...")
        end

        for t in 1:tmax # Iterate over time steps

            if verbose && (t-1)%ceil(Int, tmax/4)==0
                println("\t"^(v_lvl+2)*"Reading time t=$t out of $tmax")
            end

            for j in 1:jmax
                for i in 1:imax
                    # NOTE: Shouldn't I be iterating over in the inner loop?
                    field[i, j, t, fieldi] = parse(readline(f))
                end
            end
        end
    end

    close(f)

    # Add field to grid
    if verbose
        println("\t"^v_lvl*"Generating vtk files...")
    end

    for t in 1:tmax # Iterates over time steps

        if verbose && (t-1)%ceil(Int, tmax/10)==0
            println("\t"^(v_lvl+2)*"Saving vtk time t=$t out of $tmax")
        end

        for fieldi in 1:fieldmax # Iterates over fields

            field_name = header[fieldi]
            field_type = "scalar"
            entry_type = "node"
            field_data = reshape(field[:, :, t, fieldi], grid.nnodes)

            gt.add_field(grid, field_name, field_type, field_data,
                            entry_type; raise_warn=false)
        end

        # Save vtk of this time step
        gt.save(grid, outputname; path=save_path, num=t-1)

    end

    if paraview
        str = joinpath(save_path, outputname*"...vtk;")
        run(`paraview --data=$str`)
    end
end
