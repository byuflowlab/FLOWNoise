#=##############################################################################
# DESCRIPTION
    Converting time-pressure history to sound pressure level vs frequency and OASPL

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################


"""
    `pressure_time2frequency(time, pressure, windowtype, SPLtype)`

Converts pressure from time domain to frequency domain, calculates overall sound pressure level (OASPL)

# ARGUMENTS
* `time::Array{Real, 1}`        : Vector of time segments associated with each pressure.
* `pressure::Array{Real, 1}`    : Pressure as a function of time.
                                    
# OPTIONAL ARGUMENTS
* `windowtype::String`          : Data windowing options:
                                    "None", "Hann"
* `windowtype::String`          : Specifies wether non-weighted or A-weighted SPL should be calculated. Options:
                                    "Normal", "A-weighted"

"""
function pressure_time2frequency(time::Array{Real, 1}, pressure::Array{Real, 1}, 
                                    windowtype::String="None",
                                    SPLtype::String="A-weighted")
    
    #Define constants
    K1 = 2.243*10^16
    K3 = 1.562
    f1 = 20.599
    f2 = 107.653
    f3 = 737.862
    f4 = 12194.22 
    
    #Time & frequency parameters
    deltat = time[2] - time[1] #time spacing
    T = time[end] #total time
    deltaf = 1/T #frequency spacing
    N = length(time) # number of time points
    M = convert(Int64,floor(N/2 + 1)) # number of frequency bins of type Int64 for indexing later

    # Window data if specified - treats data as periodic for fourier transform
    if windowtype != "None"
        pressure = windowdata!(time, pressure, windowtype)
    end

    # Use FFTW package to perform fast fourier transform on pressure data
    Y = FFTW.fft(pressure)

    # Initialize variables
    Pc = zeros(Complex{Float64},M,1)
    p2local = zeros(Float64,M,1)
    SPLlocal_dB = zeros(Float64,M,1)
    p2Alocal = zeros(Float64,M,1)
    SPLlocal_dBA = zeros(Float64,M,1)
    wc = ones(Float64,M,1)
    wa = ones(Float64,M,1)
    p2 = 0.0

    x = 1:M
    f = x*deltaf

    for m in 1:M

        if m > 1
            Pc[m] = Y[m] * 2 * deltat / T #Complex pressure (magnitude and phase) at each frequency bin
        else
            Pc[m] = Y[m] * deltat / T
        end

        if SPLtype == "A-weighted" #if "Normal" then wa[m] = 1.0 and does not affect output

            #A-weighted prep, from PSU-WOPWOP manual
            wc[m] = K1 * f[m]^4 / ((f[m]^2 + f1^2)^2*(f[m]^2 + f4^2)^2)
            wa[m] = wc[m] * K3 * f[m]^4 / ((f[m]^2 + f2^2)*(f[m]^2 + f3^2))

        elseif SPLtype != "Normal" 
            error("Invalid SPLtype specified.")
        end

        # global p2 += wa[m] * abs2(Pc[m]) #overall mean square pressure
        p2 += wa[m] * abs2(Pc[m]) #overall mean square pressure
        p2local[m] = 0.5 * wa[m] * abs2(Pc[m]) #local mean square pressure for specific frequency
        SPLlocal_dB[m] = 10 * log(10,p2local[m]/((20*10^-6)^2)) #local SPL

    end

    p2 *= 0.5

    OASPL_dB = 10 * log(10,p2/(20*10^-6)^2) # overall sound pressure level in dB 

    return SPLlocal_dB, OASPL_dB
end


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
function windowdata!(time::Array{Real, 1}, data::Array{Real, 1}, windowtype::String="Hann")
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

#TODO: currently not solving problem in case 2 not matching up
end