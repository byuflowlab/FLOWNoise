#=##############################################################################
# DESCRIPTION
    Time-domain processing tools

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Sept 2020
  * License   : MIT
=###############################################################################



"""
    Read MATLAB binary file
"""
function loadbin(path::String, IDname::String, IDnum::Int, CHnum::Int;
                                            N=Inf, Nstart=0, readtype=Float32)

    # Open file
    filename = IDname*@sprintf("%03.0f", IDnum)*"_"*@sprintf("%03.0f", CHnum)*".bin"
    f = open(joinpath(path, filename), "r")

    # Determine file size
    file_size = stat(f).size
    nentries = file_size/sizeof(readtype)

    if nentries%1 == 0
        nentries = Int(nentries)
    else
        error("Invalid type for this file size: Number of entries $(nentries) is not integer")
    end

    # Move to requested start
    read(f, readtype, Nstart)

    # Determine number of entries to read
    Nread = N==Inf ? nentries - Nstart : N

    # Read data
    out = read(f, readtype, Nread)

    close(f)
    return out
end

"""
    `crosscorrelation(ps1, ps2, shift::Int, dt)`

Correlation between signal ps1 and ps2 evaluated at t=dt*shift.
"""
function crosscorrelation(ps1, ps2, shift::Int, dt)
    n1, n2 = length(ps1), length(ps2)

    if shift<0
        mini = abs(shift)+1
        maxi = min(n1, n2)
    else
        mini = 1
        maxi = min(n1, n2-shift)
    end

    val = 0
    for i in mini:maxi
        val += ps1[i]*ps2[i+shift]
    end

    return dt*val
end

"""
    `convolution(h, ps, ti::Int, dtau)`

Convolution between transfer function `h` and signal `ps` evaluated at
tau=dtau*ti.
"""
function convolution(h, ps, ti::Int, dtau)

    mini = (ti+1)-length(ps)
    mini = mini<1 ? 1 : mini
    maxi = ti

    val = 0
    for taui in mini:maxi
        val += h[taui]*ps[(ti+1)-taui]
    end

    return dtau*val
end
