#=##############################################################################
# DESCRIPTION
    Post-processing acoustic methods.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################


"Adds two sound pressure levels together."
function addSPL(x, y)

    #* For better understanding of what's going on, the following commented
    #* steps are equivalent to the returned line
    # pref = 2e-5

    # rms1 = pref^2 * 10 ^ (oaspl1 / 10) # convert to root mean square pressure
    # rms2 = pref^2 * 10 ^ (oaspl2 / 10)
    # rms_combined = rms1 + rms2 # add together

    # oaspl_combined = 10 * log10(rms_combined / (pref^2)) # convert back to decibels

    # return oaspl_combined

    return 10 * log10.(10 .^ (x/10) .+ 10 .^ (y/10))
end

"Adds a vector of sound pressure levels."
function addSPL(x::AbstractArray)

    return 10 * log10.(sum(10 .^ (x/10)))
end


"SPL spectrum to OASPL - works for A-weighted and non-A-weighted
    spectrums"
function SPL2OASPL(spls)                # vector of sound pressure levels over the spectrum

    # don't need a reference pressure when converting back and forth
    ps = 10 .^ (spls / 10)
    ps_sum = sum(ps)

    oaspl = 10 * log10(ps_sum)

    return oaspl
end


"Takes a sound pressure level spectrum in the frequency domain and
    converts it to pressure in the time domain through an inverse
    fast fourier transform."
function spl2pressure(spls)             # vector of sound pressure levels over the spectrum

    pref = 20e-6 # reference pressure
    p2s = pref^2 * 10 .^ (spls / 10.0)

    # inverse fourier transform
    ps = FFTW.ifft(p2s)

    return ps
end


"Takes in an unweighted SPL in decibels and frequency and returns
    the A-weighted SPL in decibels"
function aWeight(freq,                  # frequency
                spl)                    # unweighted sound pressure level at each frequency

    if size(freq, 1) != size(spl, 1) # ensure vectors have same dimensions
        spl = spl'
    end

    # unchanging constants defined in wopwop user manual
    K1 = 2.243e16 #s^-4
    K3 = 1.562
    f1 = 20.599 # Hz
    f2 = 107.653 # Hz
    f3 = 737.862 # Hz
    f4 = 12194.22 # Hz

    pref = 2e-5 #? may not need?

    wc = K1 * freq .^ 4 ./ ((freq .^ 2 + f1^2) .^ 2 .* (freq .^ 2 + f4^2) .^ 2) # C-weighting function
    wa = wc .* K3 .* freq .^ 4 ./ ((freq .^ 2 + f2^2) .* (freq .^ 2 + f3^2)) # A-weighting function

    p = pref^2 * 10.0 .^ (spl / 10) # convert to Pa
    pa = wa .* p # a-weight

    spla = 10 * log10.(pa / pref^2) # back to dB

    return spla
end

"Takes pressure output from wopwop and converts it to sound pressure
    level spectrum or a scalar overall sound pressure level if
    oasplFlag == true (for both cases, it is A-weighted if
    AweightingFlag == true).

Note 2: Currently does not include any windowing or filtering methods.
"
function pressure2SPL(time,                     # time history
                        p;                      # pressure (Pa) at each time step
                        AweightingFlag=false,   # SPL values are A-weighted
                        oasplFlag=false)        # returns scalar OASPL en lieu of SPL spectrum

    pref = 2e-5

    deltat = time[2] - time[1]
    T = time[end]
    deltaf = 1/T
    N = length(time)
    M = convert(Int64, floor(N/2 + 1))

    x = 0:M-1
    f = x * deltaf # frequencies

    Y = FFTW.fft(p)
    Y = Y[1:M] # cutoff frequencies higher than Nyquist frequency

    Pc = Y * 2 * deltat / T
    Pc[1] *= 0.5 # complex pressure

    p2local = 0.5 * abs2.(Pc) #mean square pressure
    spls = 10 * log10.(p2local / pref^2) #sound pressure level

    if AweightingFlag == true
        spls = aWeight(f, spls)
        if oasplFlag
            p2local = pref^2 * 10 .^ (spls / 10)
        end
    end

    if oasplFlag
        p2 = 1.0 * sum(p2local)
        oaspl = 10 * log10(p2 / pref^2)

        return oaspl
    else
        return f, spls
    end
end


"Same as function above, but weights the SPLs based on the frequency bin width.

Note: Intuitively this seems like it should work better. It will not match OASPL output from PSU-WOPWOP."
function pressure2SPL_weighted(time,                     # time history
                                p;                      # pressure (Pa) at each time step
                                AweightingFlag=false,   # SPL values are A-weighted
                                oasplFlag=false)        # returns scalar OASPL en lieu of SPL spectrum

    pref = 2e-5

    deltat = time[2] - time[1]
    T = time[end]
    deltaf = 1/T
    N = length(time)
    M = convert(Int64, floor(N/2 + 1))

    x = 0:M-1
    f = x * deltaf # frequencies

    Y = FFTW.fft(p)
    Y = Y[1:M] # cutoff frequencies higher than Nyquist frequency

    Pc = Y * 2 * deltat / T
    Pc[1] *= 0.5 # complex pressure

    p2local = 0.5 * abs2.(Pc) * deltaf #mean square pressure
    spls = 10 * log10.(p2local / pref^2) #sound pressure level

    if AweightingFlag == true
        spls = aWeight(f, spls)
        if oasplFlag
            p2local = pref^2 * 10 .^ (spls / 10)
        end
    end

    if oasplFlag
        p2 = 1.0 * sum(p2local)
        oaspl = 10 * log10(p2 / pref^2)

        return oaspl
    else
        return f, spls
    end
end


# """
#     `pressure_time2frequency(time, pressure, windowtype, SPLtype)`

# Converts pressure from time domain to frequency domain, calculates overall sound pressure level (OASPL)

# # ARGUMENTS
# * `time::Array{Real, 1}`        : Vector of time segments associated with each pressure.
# * `pressure::Array{Real, 1}`    : Pressure as a function of time.

# # OPTIONAL ARGUMENTS
# * `windowtype::String`          : Data windowing options:
#                                     "None", "Hann"
# * `windowtype::String`          : Specifies wether non-weighted or A-weighted SPL should be calculated. Options:
#                                     "Normal", "A-weighted"

# """
# function pressure_time2frequency(time::Array{Float64, 1}, pressure::Array{Float64, 1},
#                                     windowtype::String, SPLtype::String)

#     #Define constants
#     K1 = 2.243*10^16
#     K3 = 1.562
#     f1 = 20.599
#     f2 = 107.653
#     f3 = 737.862
#     f4 = 12194.22

#     #Time & frequency parameters
#     deltat = time[2] - time[1] #time spacing
#     T = time[end] #total time
#     deltaf = 1/T #frequency spacing
#     N = length(time) # number of time points
#     M = convert(Int64,floor(N/2 + 1)) # number of frequency bins of type Int64 for indexing later

#     # println("Pressure before windowing:")
#     # println(pressure)
#     # Window data if specified - treats data as periodic for fourier transform
#     if windowtype != "None"
#         pressure = FLOWNoise.windowdata!(time, pressure, windowtype)
#     end

#     # println("Pressure after windowing:")
#     # println(pressure)

#     # Use FFTW package to perform fast fourier transform on pressure data
#     Y = FFTW.fft(pressure)

#     # Initialize variables
#     Pc = zeros(Complex{Float64},M,1)
#     p2local = zeros(Float64,M,1)
#     SPLlocal_dB = zeros(Float64,M,1)
#     p2Alocal = zeros(Float64,M,1)
#     SPLlocal_dBA = zeros(Float64,M,1)
#     wc = ones(Float64,M,1)
#     wa = ones(Float64,M,1)
#     p2 = 0.0

#     x = 1:M
#     f = x*deltaf

#     for m in 1:M

#         if m > 1
#             Pc[m] = Y[m] * 2 * deltat / T #Complex pressure (magnitude and phase) at each frequency bin
#         else
#             Pc[m] = Y[m] * deltat / T
#         end

#         if SPLtype == "A-weighted" #if "Normal" then wa[m] = 1.0 and does not affect output

#             #A-weighted prep, from PSU-WOPWOP manual
#             wc[m] = K1 * f[m]^4 / ((f[m]^2 + f1^2)^2*(f[m]^2 + f4^2)^2)
#             wa[m] = wc[m] * K3 * f[m]^4 / ((f[m]^2 + f2^2)*(f[m]^2 + f3^2))

#         elseif SPLtype != "Normal"
#             error("Invalid SPLtype specified.")
#         end

#         # global p2 += wa[m] * abs2(Pc[m]) #overall mean square pressure
#         p2 += wa[m] * abs2(Pc[m]) #overall mean square pressure
#         p2local[m] = 0.5 * wa[m] * abs2(Pc[m]) #local mean square pressure for specific frequency
#         SPLlocal_dB[m] = 10.0 * log(10,p2local[m]/((20*10.0^-6)^2)) #local SPL

#     end

#     p2 *= 0.5

#     OASPL_dB = 10.0 * log(10,p2/(20*10.0^-6)^2) # overall sound pressure level in dB

#     return SPLlocal_dB, OASPL_dB
# end


"""
    `windowdata!(time, data, windowtype)`

Windows data, or "tapers" the ends to zero. This effectively makes the data periodic to avoid discontinuities in a fourier transform.

# ARGUMENTS
* `time::Array{Real, 1}`        : Vector of time segments associated with each data point in data.
* `data::Array{Real, 1}`        : Data (later to be transformed) as a function of time.

# OPTIONAL ARGUMENTS
* `windowtype::String`          : Data windowing options:
                                    "Hann", "Flat Top", "Blackman"

NOTES:
* Windowing changes the amplitude of the data in the time history, so a scaling factor must also be used, depending on which type of windowing is being performed.
* Windowing adds certain frequencies to the spectrum. Do not trust windoewed time signals which have frequencies closer than five bins away from both ends of the spectrum's frequency limits.
"""
function windowdata!(time, data, windowtype::String="Hann")
    T = time[end]
    tindex = 1
    # Apply windowing
    for t in time
        if windowtype == "Hann"
            wh = 1/2 - 1/2 * cos(2*pi*t/T)
            scalingfactor = 1.633
        elseif windowtype == "Flat Top"
            wh = 1 - 1.93 * cos(2*pi*t/T) + 1.29 * cos(4*pi*t/T) - 0.388 * cos(6*pi*t/T) + 0.0322 * cos(8*pi*t/T)
            scalingfactor = 0.515
        elseif windowtype == "Blackman"
            wh = 0.42 - 0.50 * cos(2*pi*t/T) + 0.08 * cos(4*pi*t/T)
            scalingfactor = 1.8119
        else
            error("Invalid windowtype defined.")
        end

        data[tindex] *= wh * scalingfactor
        # global tindex += 1
        tindex += 1
    end # TODO: broadcast this in the future

    return data

#TODO: currently not solving problem in case 2 not matching up
end
