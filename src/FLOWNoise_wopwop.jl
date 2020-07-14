#=##############################################################################
# DESCRIPTION
    Tools for formatting PSU-WOPWOP's inputs/outputs.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################


"""
    `read_wopwopoutput(outputname; read_path="", tec=false)`

Reads an output PSU-WOPWOP field of observer data. It returns `(header, field)`
where `header[k]` is a string with the name of the k-th field, and
`field[i, j, t, k]` is the value at the node of coordinates (i, j) at time t
in the k-th field.

Give it `tec=true` to read a TecPlot file (output of a single observer), and it
will return `(header, field)` where `field[t, k]` is the value at time t in the
k-th field.

# ARGUMENTS
* `outputname::String`           : PSU-WOPWOP output file to read. For instance,
                                    `outputname=pressure` will look for the files
                                    pressure.nam and pressure.fn under path
                                    `read_path`.

# OPTIONAL ARGUMENTS
* `read_path::String`           : Path from where to read the inputs files.
"""
function read_wopwopoutput(outputname::String; read_path="", verbose=true,
                            v_lvl=0, tec=false)

    if tec
        f = open(joinpath(read_path, outputname*".tec"), "r")

        readline(f) # Title
        header = split(replace(readline(f)[13:end], "\"" => ""), ", ")
        readline(f) # Info
        readline(f) # Empty line

        # Data
        data = []

        for ln in eachline(f)

            # Get rid of space at start of line
            while ln[1] == ' '
                ln = ln[2:end]
            end

            # Get rid of space at end of line
            while ln[end] == ' '
                ln = ln[1:end-1]
            end

            # Get rid of double spaces
            while occursin("  ", ln)
                ln = replace(ln, "  " => " ")
            end

            # Parse data
            push!(data, Meta.parse.(split(ln, " ")))
        end

        field = hcat(data...)'

        close(f)
    else
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
                header[i] = replace(header[i], str => "")
            end
        end

        # Read function file
        f = open(joinpath(read_path, outputname*".fn"), "r")
        imax, jmax, tmax, fieldmax = [Meta.parse(elem) for elem in split(readline(f))]

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
                        field[i, j, t, fieldi] = Meta.parse(readline(f))
                    end
                end
            end
        end

        close(f)
    end

    return header, field
end

"""
    `read_wopwoploading(filename; read_path="", verbose=false, v_lvl=0)`

Read a PSU-WOPWOP functional file containing loading data (v1.0). Returns
`(hdr, names, timeinfo, dims, time, data)` where data[k][:, i, j] is the value at
the i-th node (or face) at the j-th time in the k-th zone
"""
function read_wopwoploading(filename; read_path="", verbose=false, v_lvl=0)

    # Open PSU-WOPWOP file
    f = open(joinpath(read_path, filename), "r")

    # Read format header
    header = []
    push!(header, read(f, Int32, 1))                #1 Magic number
    push!(header, read(f, Int32, 2))                #2 Version
    chars = [Char(c) for c in read(f, Int8, 1024)]  #3 Comments
    str = ""; for c in chars; str *= c; end;
    push!(header, str)
    push!(header, read(f, Int32, 1))                #4 Functional file
    nz = read(f, Int32, 1)[1]                       #5 Number of zones
    push!(header, nz)
    push!(header, read(f, Int32, 1))                #6 1==structured, 2==unstructured
    push!(header, read(f, Int32, 1))                #7 1==Constant, 2==Periodic, 3==Aperiodic
    ncnt = read(f, Int32, 1)[1]                     #8 1==node-centered, 2==face-centered
    push!(header, ncnt)
    push!(header, read(f, Int32, 1))                #9 1==pressure, 2==loading, 3==flow
    push!(header, read(f, Int32, 1))                #10 1==stationary ground fixed, 2==rotating ground fixed, 3==patch fixed
    push!(header, read(f, Int32, 1))                #11 1==single, 2==double
    push!(header, read(f, Int32, 1))                #12 something
    push!(header, read(f, Int32, 1))                #13 something else

    if verbose
        println("\t"^(v_lvl)*"$(header[1])\t # Magic number")
        println("\t"^(v_lvl)*"$(header[2])\t # Version")
        println("\t"^(v_lvl)*"$(header[3])\t # Comment")
        println("\t"^(v_lvl)*"$(header[4])\t # Functional file")
        println("\t"^(v_lvl)*"$(header[5])\t\t # Number of zones")
        println("\t"^(v_lvl)*"$(header[6])\t # 1==structured, 2==unstructured")
        println("\t"^(v_lvl)*"$(header[7])\t # 1==Constant, 2==Periodic, 3==Aperiodic")
        println("\t"^(v_lvl)*"$(header[8])\t\t # 1==node-centered, 2==face-centered")
        println("\t"^(v_lvl)*"$(header[9])\t # 1==pressure, 2==loading, 3==flow")
        println("\t"^(v_lvl)*"$(header[10])\t # 1==stationary ground, 2==rotating ground, 3==patch")
        println("\t"^(v_lvl)*"$(header[11])\t # 1==single, 2==double")
        println("\t"^(v_lvl)*"$(header[12])\t # something")
        println("\t"^(v_lvl)*"$(header[13])\t # something else")
    end

    # Error cases: Current implementation only takes unstructured, constant,
    # face centered patches
    Int(header[2][1]) != 1 ? error("Only v1.0 is supported") :
    Int(header[4][1]) != 2 ? error("File is not flagged as a functional file") :
    # Int(header[6][1]) != 2 ? error("Only unstructured patches are supported") :
    # Int(header[7][1]) != 1 ? error("Only constant patches are supported") :
    # Int(header[8][1]) != 2 ? error("Only face-centered patches are supported") :
    Int(header[9][1]) != 2 ? error("Only loading vectors are supported") :
    Int(header[11][1]) != 1 ? error("Only single-precision floats are supported") :
    nothing;

    # Zone specification
    zone_specs = zeros(Int, read(f, Int32))     # Read number of zones
    for i in 1:length(zone_specs)               # Read each zone number
        zone_specs[i] = Int(read(f, Int32))
    end

    if verbose
        println("\t"^(v_lvl+1)*"Zone specification: $zone_specs")
    end

    names = String[]
    timeinfo = []
    dims = []
    ndata = zeros(Int, length(zone_specs))
    nt = nothing

    # Data header
    for zi in 1:length(zone_specs)

        # Name
        chars = [Char(c) for c in read(f, Int8, 32)]
        str = ""; for c in chars; str *= c; end;
        push!(names, str)

        # Time information
        if Int(header[7][1])==2         # Periodic case
            push!(timeinfo, Any[read(f, Float32), read(f, Int32)])
            ntimes = timeinfo[end][2]

        elseif Int(header[7][1])==3     # Aperiodic case
            push!(timeinfo, read(f, Int32))
            ntimes = timeinfo[end]

        else                            # Constant case
            push!(timeinfo, 1)
            ntimes = 1
        end

        if nt != nothing && nt != ntimes
            error("Logic error: Got different times ($(nt)!=$(ntimes))")
        end
        nt = ntimes

        # Dimensions
        if Int(header[6][1])==1         # Structured case
            push!(dims, read(f, Int32, 2))
            ndata[zi] = dims[end][1]*dims[end][2]
        else
            push!(dims, read(f, Int32))
            ndata[zi] = dims[end]
        end

        if verbose
            println("\t"^(v_lvl+1)*"$(replace(names[end], " " => ""))"*
                    "\t$(timeinfo[end])\t$(dims[end]) # Name, Time info, dims")
        end
    end

    # Time of every data block of every zone
    time = [zeros(nt) for i in 1:length(zone_specs)]

    # Data of every zone, where data[k][:, i, j] is the value at the i-th node
    # (or face) at the j-th time in the k-th zone
    data = [zeros(3, nd, nt) for nd in ndata]

    for ti in 1:nt                          # Iterate over time entries
        for zi in 1:length(zone_specs)      # Iterate over zones

            if Int(header[7][1])!=1
                time[zi][ti] = Float64(read(f, Float32))
            end

            # NOTE: Here I assume is loading data (three-dimensional)
            for dim in 1:3                  # Iterate over spatial dimension
                for di in 1:ndata[zi]       # Iterate over nodes / faces

                    data[zi][dim, di, ti] = Float64(read(f, Float32))

                end
            end
        end
    end


    hdr = Dict(
                "version"   => header[2],
                "comment"   => header[3],
                "zones"     => Int(header[5][1]),
                "structured"=> Int(header[6][1])==1,
                "time"      => Int(header[7][1]),
                "ntimes"    => nt,
                "node-centered"=> Int(header[8][1])==1,
                "data-type" => Int(header[9][1]),
                "frame"     => Int(header[10][1])
              )

    return hdr, names, timeinfo, dims, time, data
end


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
                        pref="",
                        read_path="", prompt=true,
                        verbose=true, v_lvl=0, createpath=true,
                        paraview=false)

    if createpath
        gt.create_path(save_path, prompt)
    end

    header, field = read_wopwopoutput(outputname; read_path=read_path,
                                                verbose=verbose, v_lvl=v_lvl)

    tmax = size(field, 3)
    fieldmax = size(field, 4)

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
        gt.save(grid, pref*outputname; path=save_path, num=t-1)

    end

    str = pref*outputname*"...vtk;"

    if paraview
        run(`paraview --data=$(save_path)/$(str)`)
    end

    return str
end

"""
    `geomwopwop2vtk(filename; read_path="", save_path=nothing,
pref_save_fname::String="automatic", paraview=false, verbose=false, v_lvl=0)`

Converts an input patch file in PSU-WOPWOP format to VTK format.
"""
function geomwopwop2vtk(filename; read_path="", save_path=nothing,
                        pref_save_fname::String="automatic", paraview=false,
                        verbose=false, v_lvl=0, loading_file=nothing)

    # String carrying vtk files
    if save_path!=nothing && paraview
        vtk_str = save_path*"/"
    end

    # Build output file name prefix
    if pref_save_fname=="automatic"
        doti = findfirst(x->x=='.', filename)
        if length(doti)!=0
            # Get rid of extension
            preff = filename[1:doti[end]-1]
        else
            preff = filename
        end
    else
        preff = pref_save_fname
    end

    # Open loading file
    if loading_file != nothing
        (load_hdr, load_names,
        load_timeinfo, load_dims,
        load_time,
        load_data) = read_wopwoploading(loading_file; read_path=read_path,
                                            verbose=verbose, v_lvl=v_lvl+1)
    end

    # Open PSU-WOPWOP file
    f = open(joinpath(read_path, filename), "r")

    # Read format header
    header = []
    push!(header, read(f, Int32, 1))                #1 Magic number
    push!(header, read(f, Int32, 2))                #2 Version
    chars = [Char(c) for c in read(f, Int8, 32)]    #3 Units
    str = ""; for c in chars; str *= c; end;
    push!(header, str)
    chars = [Char(c) for c in read(f, Int8, 1024)]  #4 Comments
    str = ""; for c in chars; str *= c; end;
    push!(header, str)
    push!(header, read(f, Int32, 1))                #5 Geometry file
    nz = read(f, Int32, 1)[1]                       #6 Number of zones
    push!(header, nz)
    push!(header, read(f, Int32, 1))                #7 1==structured, 2==unstructured
    push!(header, read(f, Int32, 1))                #8 1==Constant, 2==Periodic, 3==Aperiodic
    ncnt = read(f, Int32, 1)[1]                     #9 Normal centered 1==node, 2==face
    push!(header, ncnt)
    push!(header, read(f, Int32, 1))                #10 1==single, 2==double
    push!(header, read(f, Int32, 1))                #11 iblank
    push!(header, read(f, Int32, 1))                #12 something

    if verbose
        println("\t"^(v_lvl)*"$(header[1])\t # Magic number")
        println("\t"^(v_lvl)*"$(header[2])\t # Version")
        println("\t"^(v_lvl)*"$(header[3])\t # Units")
        println("\t"^(v_lvl)*"$(header[4])\t # Comment")
        println("\t"^(v_lvl)*"$(header[5])\t # Geometry file")
        println("\t"^(v_lvl)*"$(header[6])\t\t # Number of zones")
        println("\t"^(v_lvl)*"$(header[7])\t # 1==structured, 2==unstructured")
        println("\t"^(v_lvl)*"$(header[8])\t # 1==Constant, 2==Periodic, 3==Aperiodic")
        println("\t"^(v_lvl)*"$(header[9])\t\t # Normal centered 1==node, 2==face")
        println("\t"^(v_lvl)*"$(header[10])\t # 1==single, 2==double")
        println("\t"^(v_lvl)*"$(header[11])\t # iblank")
        println("\t"^(v_lvl)*"$(header[12])\t # something")
    end

    # Error cases:
    Int(header[2][1]) != 1 ? error("Only v1.0 is supported") :
    # Int(header[7][1]) != 2 ? error("Only unstructured patches are supported") :
    # Int(header[8][1]) != 1 ? error("Only constant patches are supported") :
    # Int(header[9][1]) != 2 ? error("Only face-centered patches are supported") :
    Int(header[10][1]) != 1 ? error("Only single-precision floats are supported") :
    nothing;

    if loading_file != nothing
        if load_hdr["structured"] == (Int(header[7][1]) == 1)
            error("Geometry and loading must both be either structured or unstructured!")
        end
    end

    Tflag = Int(header[8][1])       # 1==Constant, 2==Periodic, 3==Aperiodic

    # Read patch headers
    names = []
    nbnodes = []
    nbfaces = []
    imaxs = []
    jmaxs = []
    connectivity = []
    periods = []
    nts = []
    celltypes = []

    for zonei in 1:nz
        # NOTE: Here I assume it is UNSTRUCTURED and constant

        chars = [Char(c) for c in read(f, Int8, 32)]
        str = ""; for c in chars; str *= c; end;

        push!(names, str)
        push!(connectivity, [])

        if Tflag==2
            # Period
            push!(periods, Float64(read(f, Float32)))
        end
        if Tflag==2 || Tflag==3
            # nt
            push!(nts, read(f, Int32))
        end

        # Unstructured case
        if Int(header[7][1]) == 2
            push!(nbnodes, read(f, Int32))
            push!(nbfaces, read(f, Int32))
            if verbose
                println("\t"^(v_lvl)*"$(names[end])\t$(nbnodes[end]) $(nbfaces[end])\t# Name, nbNodes, nbFaces")
            end

            for facei in 1:nbfaces[end]       # Iterate over faces (cells)

                this_nbnodes = read(f, Int32) # Number of node indices to read
                push!(connectivity[end], [])

                for nodei in 1:this_nbnodes   # Collect node indices
                    push!(connectivity[end][end], read(f, Int32))
                end
            end

            push!(celltypes, -1)

        # Structured case
        else
            push!(imaxs, read(f, Int32))
            push!(jmaxs, read(f, Int32))
            push!(nbnodes, imaxs[end]*jmaxs[end])

            if verbose
                println("\t"^(v_lvl)*"$(names[end])\t$(imaxs[end]) $(jmaxs[end])\t# Name, iMax, jMax")
            end

            # Create connectivity: Quadrilateral-faces case
            if imaxs[end]!=1 && jmaxs[end]!=1
                push!(nbfaces, (imaxs[end]-1)*(jmaxs[end]-1))
                for ni in 1:imaxs[end]-1
                    for nj in 1:jmaxs[end]-1
                        # Subscripts of every node
                        subs = ( (ni, nj), (ni+1, nj), (ni+1, nj+1), (ni, nj+1) )
                        # Convert subscripts to linear index of every node
                        this_cell = [Base._sub2ind((imaxs[end], jmaxs[end]), sub...) for sub in subs]
                        push!(connectivity[end], this_cell)
                    end
                end

                push!(celltypes, -1)

            # Create connectivity: Compact-patch (line) case
            else
                push!(nbfaces, imaxs[end]*jmaxs[end]-1)
                for ni in 1:(imaxs[end]*jmaxs[end]-1)
                    push!(connectivity[end], [ni, ni+1])
                end

                push!(celltypes, 4)
            end
        end

        if verbose && Tflag==2
            println("\t"^(v_lvl+1)*"$(periods[end])\t$(nts[end])\t# Period, nt")
        elseif verbose && Tflag==3
            println("\t"^(v_lvl+1)*"$(nts[end])\t# nt")
        end

    end


    # -------------------- READ POINTS AND NORMALS -----------------------------

    Nflag = Int(header[9][1]) # Normal centered 1==node, 2==face

    if Tflag != 1
        points = [zeros(3, nbnodes[zi], nts[zi]) for zi in 1:nz]
        normals = [zeros(3, Nflag==1 ? nbnodes[zi] : nbfaces[zi], nts[zi]) for zi in 1:nz]
        times = [zeros(nts[zi]) for zi in 1:nz]
    else
        points = [zeros(3, nbnodes[zi]) for zi in 1:nz]
        normals = [zeros(3, Nflag==1 ? nbnodes[zi] : nbfaces[zi]) for zi in 1:nz]
        times = []
    end

    # Format cells as vtks (0-indexed)
    vtk_cells = [[Int[ind-1 for ind in cell] for cell in connectivity[zi]] for zi in 1:nz]

    # Zone of every point
    point_zones = [[zi for pi in 1:nbnodes[zi]] for zi in 1:nz]

    str = ""

    # NOTE: I assume that all zones have the same number of time steps, so here
    # I'm checking that assumption
    if Tflag != 1
        for zi in 1:nz-1
            if nts[zi]!=nts[zi+1]
                error("Found that not all zones have the same number of time"*
                        " steps ($(nts[zi])!=$(nts[zi+1])).")
            end
        end
    end

    for ti in 1:(Tflag != 1 ? nts[1] : 1)   # Iterate over time steps
        for zonei in 1:nz                   # Iterate over zones

            # Read time signature
            if Tflag != 1
                times[zonei][ti] = Float64(read(f, Float32))
            end

            # Read points
            for k in 1:3
                for pi in 1:nbnodes[zonei]
                    points[zonei][k, pi, ti] = Float64(read(f, Float32))
                end
            end

            # Read normals
            for k in 1:3
                for ni in 1:(Nflag==1 ? nbnodes[zonei] : nbfaces[zonei])
                    normals[zonei][k, ni, ti] = Float64(read(f, Float32))
                end
            end

            if save_path != nothing

                # Formats points as vtks
                vtk_points = [ points[zonei][:, pi, ti] for pi in 1:nbnodes[zonei] ]

                # Format normals as a vtk field
                vtk_normals = [ normals[zonei][:, pi, ti] for pi in 1:(
                                    Nflag==1 ? nbnodes[zonei] : nbfaces[zonei]) ]

                # Data fields
                point_data = [Dict( "field_name" => "zone",
                                    "field_type" => "scalar",
                                    "field_data" => point_zones[zonei])]

                if Nflag == 2
                    cell_data = [Dict( "field_name"  => "normals",
                                        "field_type" => "vector",
                                        "field_data" => vtk_normals)]
                else
                    cell_data = nothing
                    push!(point_data, Dict( "field_name" => "normals",
                                            "field_type" => "vector",
                                            "field_data" => vtk_normals))
                end

                # Loading data
                if loading_file != nothing

                    # NOTE: Here I assume that the loading file has data for
                    # the same number of zones that the geometry file (even
                    # though that is not necessarely true, but it is for
                    # files generated in my framework)
                    loading_vectors = [load_data[zonei][:, ni, ti] for ni in 1:size(load_data[zonei], 2)]

                    loading_field = Dict(   "field_name"  => load_names[zonei],
                                            "field_type" => "vector",
                                            "field_data" => loading_vectors)

                    if load_hdr["node-centered"]
                        push!(point_data, loading_field)
                    else
                        if cell_data == nothing
                            cell_data = []
                        end
                        push!(cell_data, loading_field)
                    end
                end

                filename = preff*"_"*replace(names[zonei], " " => "")

                gt.generateVTK(filename, vtk_points; cells=vtk_cells[zonei],
                                point_data=point_data, cell_data=cell_data,
                                path=save_path, num=Tflag!=1 ? ti : nothing,
                                override_cell_type=celltypes[zonei])

                if ti==1
                    str *= filename*(Tflag==1 ? ".vtk;" : "...vtk;")
                end
            end

        end
    end

    close(f)

    if paraview && save_path != nothing
        vtk_str *= str
        run(`paraview --data=$(vtk_str)`)
    end

    return names, points, vtk_cells, periods, times, str
end

"""
Reads a time-sequency or single legacy ASCII unstructured VTK file and formats
them into a single file in WOPWOP's geometry format v1.0.
"""
function vtk2wopwop(in_filename, out_filename; read_path="", nums=nothing,
                            save_path="", wopbin=true,
                            t0=0.0, tf=1.0, period=nothing, compact=false)

    _in_filename = in_filename*(in_filename[end-3:end]==".vtk" ? "" : ".vtk")

    # Time type: 1==Constant, 2==periodic, 3==aperiodic
    nums==nothing   ? Tflag = 1 :
    period!=nothing ? Tflag = 2 :
                      Tflag = 3

    # Read first file
    (points, cells, cell_types,
        data) = gt.read_vtk(Tflag==1 ? _in_filename :
                                replace(_in_filename, ".vtk" => ".$(nums[1]).vtk");
                                path=read_path)

    nnodes = size(points, 2)
    ncells = size(cells, 1)

    # Create output file
    f = open(joinpath(save_path, out_filename), "w")

    # Binary / ASCII printing
    prnt(x) = wopbin ? write(f, x) : print(f, x)
    prntln(x) = wopbin ? write(f, x) : print(f, x, "\n")

    # Convertion to 4-bytes numbers
    # NOTE: 4 bytes = 4*8 bites = 32 bites
    fl(x) = Float32(x)
    nt(x) = Int32(x)
    # Convertion to n-bytes string
    st(x::String, n) = x * " "^(n-length(x))

    # Magic number
    prntln(nt(42))
    # Version number
    prnt(nt(1))
    prntln(nt(0))
    # Units for Tecplot
    prntln(st("N/m^2", 32))
    # Comments
    prntln(st("Geometry input file for PSU-WOPWOP (Format v1.0)\n"*
              "------------------------------------------------\n"*
              "Created by FLOWNoise (written by Eduardo Alvarez)\n"*
              "https://github.com/byuflowlab/FLOWNoise\n"*
              "Creation date: $(Dates.now())\n"*
              "Units: SI\n"*
              (
                compact ? "Format: Structured grid, node-centered" :
                          "Format: Unstructured grid, face-centered"
              ), 1024))

    # Format string
    prntln(nt(1))               # Geometry file flag
    prntln(nt(1))               # Number of zones
    prntln(nt(2^!compact))      # 1==structured, 2==unstructured
    prntln(nt(Tflag))           # Geometry 1==constant, 2==periodic, 3==aperiodic
    prntln(nt(2^!compact))      # Normal vectors 1==node, 2==face
    prntln(nt(1))               # Floating point 1==single, 2==double
    prntln(nt(0))               # iblank values 1==included, 0==not
    prntln(nt(0))               # WOPWOP secret conspiracy

    # ----------------- PATCH ----------------------------------------
    # Name
    prntln(st(compact ? "vtkcompactpatch" : "vtkpatch", 32))
    # TimeInformation
    if Tflag==2
        # Period
        prntln(fl(period))
        # nt
        prntln(nt( length(nums) ))
    elseif Tflag==3
        # nt
        prntln(nt( length(nums) ))
    end

    if compact
        # iMax
        prntln(nt( 1 ))
        # jMax
        prntln(nt( nnodes ))
    else
        # nbNodes
        prntln(nt( nnodes ))
        # nbFaces
        prntln(nt( ncells ))
        # Connectivity
        # NOTE: Here I assume that connectivity is the same in all VTK files
        for cell in cells
            prnt(nt( size(cell, 1) ))         # Number of nodes in this cell
            # for pi in reverse(cell)         # Clockwise node ordering
            for pi in cell                    # NOTE: Turns out that FLOWVLM's VTKs are already clockwise
                prnt(nt( pi+1 ))              # 1-indexed node index
            end
            if !wopbin; prntln(""); end;
        end
    end



    # ----------------- DATA OF PATCH -------------------------------

    Ns = zeros(3, compact ? nnodes : ncells)

    for ti in 1:(Tflag==1 ? 1 : length(nums))
        # Read vtk
        if ti != 1
            # NOTE: Here I assume that `nums` is given sequentally
            (points, cells, cell_types,
             data) = gt.read_vtk(replace(_in_filename, ".vtk" => ".$(nums[ti]).vtk");
                                    path=read_path)

            Ns[:] .= 0
        end

        # Calculate normals
        if compact
            # NOTE: Here I assume that the lifting line is defined sequentially
            for pi in 1:nnodes-1
                # NOTE: Here I give it the vector that connect the points
                # as the normal, crossing my fingers that PSU-WOPWOP will never
                # use the actual normal direction for anything (I think that
                # it only expects the magnitude to be the length of this segment
                # of the line)
                D = points[:, pi+1] - points[:, pi]

                # Here I split the length in two for each node
                Ns[:, pi] += D/2
                Ns[:, pi+1] += D/2
            end

            # # Here I'm normalizing the normal by the length
            # for pi in 1:nnodes
            #     Ns[:, pi] /= norm(Ns[:, pi])
            # end

        else
            for (ci, cell) in enumerate(cells)

                if length(cell) != 4
                    error("Only quadrilateral cells are supported!")
                end

                crss1 = -cross( points[:, cell[2]+1] - points[:, cell[1]+1],
                                points[:, cell[3]+1] - points[:, cell[1]+1])
                crss2 = -cross( points[:, cell[4]+1] - points[:, cell[3]+1],
                                points[:, cell[1]+1] - points[:, cell[3]+1])

                Ns[:, ci] .= crss1/2 + crss2/2
            end
        end

        # Time signature
        if Tflag != 1
            prntln(fl( (tf-t0)*(ti-1)/(length(nums)-1) ))
        end

        # nbNodes floating point x coordinates
        # nbNodes floating point y coordinates
        # nbNodes floating point z coordinates
        for k in 1:3
            for pi in 1:nnodes
                prntln(fl(points[k, pi]))
            end
        end

        if compact
            # imax × jmax floating point normal vector x coordinates
            # imax × jmax floating point normal vector y coordinates
            # imax × jmax floating point normal vector z coordinates
            for k in 1:3
                for j in 1:nnodes
                    prntln(fl(Ns[k, j]))
                end
            end

        else
            # nbFaces floating point normal vector x coordinates
            # nbFaces floating point normal vector y coordinates
            # nbFaces floating point normal vector z coordinates
            for k in 1:3
                for j in 1:ncells
                    prntln(fl(Ns[k, j]))
                end
            end
        end
    end

    close(f)

    return nothing
end



"""
Read the solution field created by PSU-WOPWOP named `fieldname` (i.e., pressure,
spl_spectrum, OASPLdB, or OASPLdBA) under `read_path`.
"""
function read_pswfield(fieldname::String, read_path; re20=false, tec=false)
    header, field = read_wopwopoutput(fieldname; read_path=read_path, tec=tec)

    if re20
        header_hash = Dict( (occursin(",re20", h) ? h[1:findfirst(",re20", h)[1]-1] : h, i)
                                                            for (i,h) in enumerate(header))
    else
        header_hash = Dict((h, i) for (i,h) in enumerate(header))
    end

    return field, header, header_hash
end


function read_pswfield(fieldnames::Array{String, 1}, read_path;
                            re20crit=["OASPLdB", "OASPLdBA"], tec=false)

    datas = Dict()

    for fieldname in fieldnames
        data = read_pswfield(fieldname, read_path; re20=fieldname in re20crit, tec=tec)

        datas[fieldname] = Dict("field"=>data[1], "hd"=>data[2], "hs"=>data[3])
    end

    return datas
end

function fetch_pswfield(read_path::String, args...;
                        # fieldnames=["pressure", "spl_spectrum", "OASPLdB", "OASPLdBA"],
                        fieldnames=["spl_spectrum", "OASPLdB", "OASPLdBA"],
                        psw_datasets=Dict(), optargs...)

    println("*"^72*"\n*\tReading dataset $read_path\n"*"*"^72)

    psw_datasets[read_path] = read_pswfield(fieldnames, read_path, args...;
                                                                    optargs...)

    return psw_datasets[read_path]
end

function fetch_pswdataset(read_path, args...; psw_datasets=Dict(), optargs...)

    if !(read_path in keys(psw_datasets))
        return fetch_pswfield(read_path, args...; psw_datasets=psw_datasets, optargs...)
    else
        return psw_datasets[read_path]
    end

end

# ------------------------- FUNCTIONS FOR BPM ----------------------------------
read_bpmoutput(args...; optargs...) = read_wopwopoutput(args...; tec=true, optargs...)
fetch_bpmdataset(args...; optargs...) = fetch_pswdataset(args...; tec=true,
                                                         fieldnames=["spl_spectrum", "splA_spectrum", "OASPLdB", "OASPLdBA"],
                                                                            optargs...)
