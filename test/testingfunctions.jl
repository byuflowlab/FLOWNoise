"""
    `sort_data(data, file_type)`

Sorts data from wopwop text file outputs.

# ARGUMENTS
* `data::Array{Real, Real}`     : Name of the file to be read, path included
* `file_type::String`           : Specifies the type of output file from WOPWOP. Options:
                                    "pressure-time"
                                    "spl_spectrum"
                                    "oaspl"
"""
function sort_data(data, file_type)

    if file_type == "pressure-time"

        time = data[:,1]
        thickness = data[:,2]
        loading = data[:,3]
        total = data[:,4]
        thicknessX = data[:,5]
        thicknessY = data[:,6]
        thicknessZ = data[:,7]
        loadingX = data[:,8]
        loadingY = data[:,9]
        loadingZ = data[:,10]
        totalX = data[:,11]
        totalY = data[:,12]
        totalZ = data[:,13]

        return time, thickness, loading, total, thicknessX, 
                thicknessY, thicknessZ, loadingX, loadingY, 
                loadingZ, totalX, totalY, totalZ

    elseif file_type == "spl-spectrum"

        frequency = data[:,1]
        thickness_dB = data[:,2]
        loading_dB = data[:,3]
        total_dB = data[:,4]
        thickness_dBA = data[:,5]
        loading_dBA = data[:,6]
        total_dBA = data[:,7]

        return frequency, thickness_dB, loading_dB, total_dB, 
                thickness_dBA, loading_dBA, total_dBA

    elseif file_type == "oaspl" || file_type == "oaspla"

        time, thickness_OASPL, loading_OASPL, total_OASPL = 
            Float64[], Float64[], Float64[], Float64[]

        time = data[1]
        thickness_OASPL = data[2]
        loading_OASPL = data[3]
        total_OASPL = data[4]

        return time, thickness_OASPL, loading_OASPL, total_OASPL

    else
        println("ERROR: Invalid file_type input")
    end
end

"""
    `prepare_test(filepath, windowtype)`

Reads all output files from wopwop test cases and calculates the same results from the post_processing 
    functions in FLOWNoise. Returns the needed variables to perform the required tests.

# ARGUMENTS
* `filepath::String`            : Name of the file to be read, path included
* `windowtype::String`          : Data windowing options:
                                    "None", "Hann"
"""
function prepare_test(filepath, windowtype)

    #Read time-pressure data
    filename = joinpath(filepath, "pressure")
    file_type = "pressure-time"
    header, data = FLOWNoise.read_wopwopoutput(filename, read_path="", verbose=true, v_lvl=0, tec=true)
    time, thickness, loading, total, 
            thicknessX, thicknessY, thicknessZ, 
            loadingX, loadingY, loadingZ, 
            totalX, totalY, totalZ = sort_data(data, file_type)

    #Read Pressure Spectrum data
    filename = joinpath(filepath, "spl_spectrum")
    file_type = "spl-spectrum"
    header, data = FLOWNoise.read_wopwopoutput(filename, read_path="", verbose=true, v_lvl=0, tec=true)
    frequency, thickness_dB_actual, loading_dB_actual, total_dB_actual, 
            thickness_dBA_actual, loading_dBA_actual, total_dBA_actual = sort_data(data, file_type)

    #Read non-weighted OASPL data
    filename = joinpath(filepath, "OASPLdB")
    file_type = "oaspl"
    header, data = FLOWNoise.read_wopwopoutput(filename, read_path="", verbose=true, v_lvl=0, tec=true)
    time_2, thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual= sort_data(data, file_type)

    #Read A-weighted OASPL data
    filename = joinpath(filepath, "OASPLdBA")
    file_type = "oaspla"
    header, data = FLOWNoise.read_wopwopoutput(filename, read_path="", verbose=true, v_lvl=0, tec=true)
    time_3, thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual = sort_data(data, file_type)

    # Convert from time to frequency domain
    SPLtype = "Normal"
    thickness_dB_calc, thickness_OASPL_calc = FLOWNoise.pressure_time2frequency(time, thickness, windowtype, SPLtype)
    loading_dB_calc, loading_OASPL_calc = FLOWNoise.pressure_time2frequency(time, loading, windowtype, SPLtype)
    total_dB_calc, total_OASPL_calc = FLOWNoise.pressure_time2frequency(time, total, windowtype, SPLtype)
    SPLtype = "A-weighted"
    thickness_dBA_calc, thickness_OASPLA_calc = FLOWNoise.pressure_time2frequency(time, thickness, windowtype, SPLtype)
    loading_dBA_calc, loading_OASPLA_calc = FLOWNoise.pressure_time2frequency(time, loading, windowtype, SPLtype)
    total_dBA_calc, total_OASPLA_calc = FLOWNoise.pressure_time2frequency(time, total, windowtype, SPLtype)

    println("Actual thickness OASPL = ", thickness_OASPL_actual)
    println("Calculated thickness OASPL = ", thickness_OASPL_calc)

    println("Actual loading OASPL = ", loading_OASPL_actual)
    println("Calculated loading OASPL = ", loading_OASPL_calc)

    println("Actual total OASPL = ", total_OASPL_actual)
    println("Calculated total OASPL = ", total_OASPL_calc)

    println("Actual thickness OASPLA = ", thickness_OASPLA_actual)
    println("Calculated thickness OASPLA = ", thickness_OASPLA_calc)

    println("Actual loading OASPLA = ", loading_OASPLA_actual)
    println("Calculated loading OASPLA = ", loading_OASPLA_calc)

    println("Actual total OASPLA = ", total_OASPLA_actual)
    println("Calculated total OASPLA = ", total_OASPLA_calc)

    println("")

    return thickness_dB_calc, loading_dB_calc, total_dB_calc,
            thickness_dBA_calc, loading_dBA_calc, total_dBA_calc,
            thickness_dB_actual, loading_dB_actual, total_dB_actual,
            thickness_dBA_actual, loading_dBA_actual, total_dBA_actual,
            thickness_OASPL_calc, loading_OASPL_calc, total_OASPL_calc, 
            thickness_OASPLA_calc, loading_OASPLA_calc, total_OASPLA_calc, 
            thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual, 
            thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual

end

"""
    `percent_error(x1, x2)`

Simple percent error calculation to clean up code.

"""
function percent_error(x1, x2)


    return abs((x1-x2)/x1)

end