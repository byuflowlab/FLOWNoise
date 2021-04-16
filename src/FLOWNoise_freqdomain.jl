#=##############################################################################
# DESCRIPTION
    Frequency-domain processing tools

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Sept 2020
  * License   : MIT
=###############################################################################


"""
    `autospec(x, fs, ns, overlap; window=hanning, density=false)`

Calculates the autospectrum of the waveform `x` using the periodigram method.

# Arguments
* `x::Array`       : (Pa) waveform.
* `fs::Real`       : (Hz) sample rate.
* `ns::Int`        : Size of blocks.
* `overlap::Real`  : Overlap of blocks between 0 and 1.
* `window::Function`   : Windowing function.
* `density::Bool`  : Returns spectrual density if true.
"""
function autospec(x, fs::Real, ns::Int, overlap::Real; window=hanning,
                                                    density=false, pref=2e-5)

    N = ns                     # Number of samples going through the FFT
    M = floor(Int, N/2+1)      # Frequencies to keep
    Gxx = zeros(M)             # Build the mean square of each frequency
    nave = 0

    delta_freq = fs/(N-1)
    freqs = [(m-1)*delta_freq for m in 1:M]

    # Iterate over blocks
    for firsti in range(1, length(x), step=ceil(Int, ns*(1-overlap)))
        endi = firsti + ns - 1

        if endi > length(x); break; end;

        # Generate window factors
        ww = window.(1:ns, ns)
        W = mean(ww.*conj(ww))  # Used to scale the ASD for energy conservation

        # Correction factor for single-sided fft (this is a copy/paste from Dr. Gee's function)
        # aux = sqrt(2*delta_freq/N/fs/W)
        aux = sqrt(2)*(1/N)*sqrt(1/W)

        # Grab the block
        this_x = x[firsti:endi]

        # Enforce zero mean
        # this_x .-= mean(this_x)

        # Window the block
        this_x = ww.*this_x

        # Corrected single-sided FFT the block
        fftx = aux*FFTW.fft(this_x)[1:M]

        # Magnitude square of FFT
#         fftx2 = (norm.(fftx)).^2
#         fftx2 = (abs.(fftx)).^2
        fftx2 = conj(fftx).*fftx

        # Build the average
        Gxx .+= fftx2
        nave += 1
    end

    Gxx /= nave

    OASPL = 20*log10(sqrt(sum(Gxx))/pref)

    # Normalize by the frequency bin to get density
    if density
        Gxx /= delta_freq
    end

    return Gxx, freqs, OASPL
end

"""
    `crossspec(x, fs, ns, overlap; window=hanning, density=false)`

Calculates the cross-spectrum of the waveforms `x` and `y` using the periodigram
method. This function is adapted from the MATLAB function by Kent Gee and Alan
Wall:

# Arguments
* `x::Array`       : (Pa) waveform.
* `fs::Real`       : (Hz) sample rate.
* `ns::Int`        : Size of blocks.
* `overlap::Real`  : Overlap of blocks between 0 and 1.
* `window::Function`   : Windowing function.
* `density::Bool`  : Returns spectrual density if true.
"""
function crossspec(x, y, fs::Real, ns::Int, overlap::Real; window=hanning,
                                                    density=true, pref=2e-5)

    # Clip data sets
    # N = Int(2^floor(log2( min(length(x), length(y)) )))
    N = min(length(x), length(y))
    x1 = x[1:N]
    x2 = y[1:N]

    # Enforce zero mean
    x1 .-= mean(x1)
    x2 .-= mean(x2)

    # Generate window factors
    ww = window.(1:ns, ns)
    W = mean(ww.*conj(ww))  # Used to scale the ASD for energy conservation

    # SPLITS DATA INTO BLOCKS
    #  Divides total data set into blocks of length ns with overlap, and
    #  windowed.
    N = ns                     # Number of samples going through the FFT
    M = floor(Int, N/2+1)      # Frequencies to keep
    Gxy = zeros(Complex, M)             # Build the mean square of each frequency
    nave = 0

    # Correction factor for single-sided fft (this is a copy/paste from Dr. Gee's function)
    scale = 2/ns/fs/W

    delta_freq = fs/(N-1)
    freqs = [(m-1)*delta_freq for m in 1:M]

    # Iterate over blocks
    for firsti in range(1, length(x1), step=ceil(Int, ns*(1-overlap)))
        endi = firsti + ns - 1

        if endi > length(x1); break; end;

        # Grab the block and window it
        this_x1 = ww.*x1[firsti:endi]
        this_x2 = ww.*x2[firsti:endi]

        # Corrected single-sided FFT the block
        fftx1 = FFTW.fft(this_x1)[1:M]
        fftx2 = FFTW.fft(this_x2)[1:M]

        # Build the average
        Gxy .+= scale*conj(fftx1).*fftx2     # Units are Pa^2/Hz
        nave += 1
    end

    Gxy /= nave

    aux = sqrt(2)*(1/N)*sqrt(1/W)
    OASPL = 20*log10(sqrt( aux^2/scale *  sum(Gxy))/pref)

    if density==false # Multiple by the frequency bin
        Gxy *= delta_freq
    end

    return Gxy, freqs, OASPL
end


"""
    `hanning(n, N)`
Hann window of period N sampled at n
"""
hanning(n::Int, N::Int) = n>N || n<1 ? 0 : sin(pi*(n-1)/(N-1))^2

triangular(n::Int, N::Int) = n>N || n<1 ? 0 : 1 - abs((n - N/2) / (N/2))

sinw(n::Int, N::Int) = n>N || n<1 ? 0 : sin(pi*n/N)

nuttall(n::Int, N::Int) = n>N || n<1 ? 0 :
                            0.355768 - 0.487396*cos(2*pi*n/N) +
                            0.144232*cos(4*pi*n/N) - 0.012604*cos(6*pi*n/N)

gaussian(n::Int, N::Int) = n>N || n<1 ? 0 : exp(-1/2*( (n-N/2) / (0.5*N/2) )^2)

expw(n::Int, N::Int) = n>N || n<1 ? 0 : exp(-abs(n-N/2)* 1/((N/2)/(60/8.69)))


"""

# Arguments
* `ps::Array{Float64, 1}`     : (Pa) pressure waveform
* `fs::Int`                   : (Hz) sample rate
"""
function spl_spectrum_from_autospectrum(ps, fs; pref=20e-6, numblocks=10,
                                                    overlap=0.5, window=hanning)

    ns = ceil(Int, length(ps)/numblocks)                   # Block size

    Gxx, freqs, OASPL = autospec(ps, fs, ns, overlap; window=window,
                                                        density=true, pref=pref)

    ww = window.(1:ns, ns)
    wwmean = mean(ww)
    df = freqs[2]-freqs[1]

    # Convert to amplitude
    ampl = sqrt.(Gxx .* df ./ wwmean^2)

    # Convert to SPL
    spls = 20*log10.(ampl./pref)

    return freqs, spls, OASPL
end


"""

# Arguments
* `ps::Array{Float64, 1}`     : (Pa) pressure waveform
* `fs::Int`                   : (Hz) sample rate
"""
function spl_spectrum_from_crossspectrum(ps, fs; numblocks=10, optargs...)

    ns = ceil(Int, length(ps)/numblocks)                   # Block size

    return spl_spectrum_from_crossspectrum(ps, fs, ns; optargs...)

end

function spl_spectrum_from_crossspectrum(ps, fs, ns; optargs...)

    return spl_crossspectrum_from_Gxy(ps, ps, fs, ns; optargs...)
end


function spl_crossspectrum_from_Gxy(ps_x, ps_y, fs, ns; pref=20e-6,
                                                    overlap=0.5, window=hanning)

    Gxy, freqs, OASPL = crossspec(ps_x, ps_y, fs, ns, overlap; window=window,
                                                        density=true, pref=pref)

    ww = window.(1:ns, ns)
    wwmean = mean(ww)
    df = freqs[2]-freqs[1]

    # Convert to amplitude
    ampl = sqrt.(Gxy .* df ./ wwmean^2)

    # Convert to SPL
    spls = 20*log10.(ampl./pref)

    return freqs, spls, OASPL
end

"""

# Arguments
* `ps::Array{Float64, 1}`     : (Pa) pressure waveform
* `fs::Int`                   : (Hz) sample rate
"""
function spl_spectrum_from_fft(ps, fs; pref=20e-6)


    dt = 1/fs                      # (s) sample period

    T = (length(ps)-1)/fs          # (s) total sampling time
    N = length(ps)                 # Number of samples

    # Fourier Transform
    ps_fft = 1/N*FFTW.fft(ps)  # FFTW doesn't multiply by 1/N for whatever reason
    ps_fftshift = FFTW.fftshift(ps_fft)

    # Define double sided frequency array
    dfreq = 1/T                 # (Hz) frequency step
    freqmax = floor(N/2)*dfreq  # (Hz) for even N, this is the Nyquist rate.
                                #      for odd N, this is 0.5*df less than the Nyquist rate

    # The DFT is N-periodic. The N+1th sample (if you went there) is the
    # same as the first sample (the Zero frequency bin).


    if mod(N, 2)==0 # Case N is even
    #=
        For even N, fmax is the nyquist rate.
        For an even N, the first sample is the 0 frequency and the (N/2
        +1)th sample is the sample at the Nyquist rate (by whatever standard that
        defines the FFT, the Nyquist rate is placed at the end of the negative
        side of the double sided array, so the positive side only goes up one bin
        less than the nyquist.
        The Nth sample is then the sample right before the 0 frequency
        (or equivalently the sampling frequency) again. because of periodicity
        the Nth sample is the same as the first negative frequency, the N-1th
        sample is the same as the second negative frequency, and so forth.
    =#
        freqs_raw = -freqmax:dfreq:(freqmax-dfreq)

    else # Case N is odd
    #=
        For odd N fmax is not the Nyquist, its half a bin below the Nyquist
        for odd N the first sample is the zero frequency and the floor(N/2)th
        sample represents the frequcny half a bin below the Nyquist. The
        floor(N/2)+1 to Nth frequencies then, by periodicity, are the same as the
        -floor(N/2)th ascending to the sample right before the zero frequency
        again.
    =#
        freqs_raw = -freqmax:dfreq:freqmax

    end

    freqs_unshifted = (0:(N-1))*dfreq
    freqs_singlesided = 0:dfreq:freqmax


    # (multiply the complex values by sqrt(2), so that when it's
    #  squared it's the same as being multiplied by 2, because we
    #  lost half the energy by only going single sided)
    ps_fft_singlesided = sqrt(2)*ps_fft[1:(floor(Int, N/2)+1)];

    spls = 20*log10.(abs.(ps_fft_singlesided)/pref)

    OASPL_freq_singlesided = 10*log10(sum(abs.(ps_fft_singlesided).^2)./pref^2)

    return freqs_singlesided, spls, OASPL_freq_singlesided
end

# """
# Convert the spectral density Gxx from narrow band to one-third octave band.
# """
# function narrow2onethird_Gxx(freqs, Gxx)
#
#     df = freqs[2]-freqs[1]                                # (Hz) Bin width
#
#     fc = [1, 1.25, 1.6, 2.0, 2.5, 3.15, 4, 5, 6.3, 8]     # (Hz) band center frequencies
#     fc = vcat(fc, fc*10, fc*100, fc*1000, fc*10000, 1e5)
#
#     flow = [0.88, 1.13, 1.414, 1.76, 2.25, 2.825, 3.53, 4.4, 5.65, 7.07] # (Hz) lower bound of 1/3-oct bins
#     flow = vcat(flow, flow*10, flow*100, flow*1000, flow*10000, 0.88e5)
#
#     fhi = [1.13, 1.414, 1.76, 2.25, 2.825, 3.53, 4.4, 5.65, 7.07, 8.8] # (Hz) upper bound of 1/3-oct bins
#     fhi = vcat(fhi, fhi*10, fhi*100, fhi*1000, fhi*10000, 1.13e5)
#
#     freqs_oto = fc                                        # (Hz) One-third octave frequency bins
#     Gxx_oto = zeros(eltype(Gxx), length(freqs_oto))       # One-third octave autospectrum
#
#     for i in 1:length(freqs)
#
#         # Find oto bin
#         otoi = nothing
#
#         for j in 1:length(freqs_oto)
#             if flow[j] <= freqs[i] &&  freqs[i] < fhi[j]
#                 otoi = j
#                 break
#             end
#         end
#
#         if otoi != nothing
#             Gxx_oto[otoi] += Gxx[i]*df
#         end
#     end
#
#     return freqs_oto, Gxx_oto
# end

# """
# # Convert spectral SPL from narrow band to one-third octave band.
# # """
# function narrow2onethird_spl_DEPRE(freqs, spls; pref=20e-6)
#
#     # Convert SPL to amplitude
#     ampl = pref * 10 .^ (spls ./ 20)
#
#     # Convert amplitude to pseudo-density
#     df = freqs[2]-freqs[1]
#     Gxx = ampl.^2 ./ df
#
#     # Convert narrow band to one-third octave band
#     freqs_oto, Gxx_oto = narrow2onethird_Gxx(freqs, Gxx)
#
#     # 1/3 octave band bounds
#     fc = freqs_oto
#     flow = vcat(fc[2:end]/sqrt(2), fc[end-2]*sqrt(2))
#     fhi = vcat(fc[3]/sqrt(2), fc[1:end-1]*sqrt(2))
#
#     # Convert pseudo-density to amplitude
#     # dfs = vcat(freqs_oto[2]-freqs_oto[1], [freqs_oto[i+1] - freqs_oto[i] for i in 1:length(freqs_oto)-1])
#     # ampl_oto = sqrt.(Gxx_oto .* dfs)
#     # dfs = [freqs_oto[i+1] - freqs_oto[i] for i in 1:length(freqs_oto)-1]
#     # ampl_oto = sqrt.(Gxx_oto[1:end-1] .* dfs)
#
#     dfs = fhi .- flow
#     ampl_oto = sqrt.(Gxx_oto .* dfs)
#
#     # Convert amplitude to SPL
#     spls_oto = 20*log10.(ampl_oto./pref)
#
#     # return freqs_oto[1:end-1], spls_oto
#     return freqs_oto, spls_oto
# end


"""
Convert spectral SPL from narrow band to one-third octave band.
"""
function narrow2onethird_spl(freqs, spls)

    spls_oto = zeros(eltype(spls), length(fc_oto)) # One-third octave autospectrum

    # Convert SPLs to pressure
    # and distribute all narrow band bins onto OTO bins
    otoi = 1
    for (i, freq) in enumerate(freqs) # Iterate over narrow bins

        # Find oto bin
        if freq < flow_oto[1]
            otoi = 1
        elseif freq >= fhi_oto[end]
            otoi = length(fc_oto)
        else
            for j in otoi:length(fc_oto)
                if flow_oto[j] <= freq &&  freq < fhi_oto[j]
                    otoi = j
                    break
                end
                if j==length(fc_oto)
                    error("LOGIC ERROR: OTO bin not found!")
                end
            end
        end

        # println("$i\t$otoi\t$(freq)\t$(fc_oto[otoi])")

        # Convert SPL to pressure and add to this OTO bin
        spls_oto[otoi] += 10^(spls[i]/10)
    end

    # Convert OTO pressure to OTO SPL
    for i in 1:length(spls_oto)
        if spls_oto[i] == 0
            spls_oto[i] = -55
        else
            spls_oto[i] = 10*log10(spls_oto[i])
        end
    end

    return deepcopy(fc_oto), spls_oto
end


# """
# Convert spectral SPL from narrow band to one-third octave band.
# """
# function narrow2onethird_spl_DEPRE(freqs, spls; pref=20e-6)
#
#     # Convert SPL to amplitude
#     ampl = pref * 10 .^ (spls ./ 20)
#
#     # Convert amplitude to pseudo-density
#     df = freqs[2]-freqs[1]
#     Gxx = ampl.^2 ./ df
#
#     # Convert narrow band to one-third octave band
#     freqs_oto, Gxx_oto = narrow2onethird_Gxx(freqs, Gxx)
#
#     # 1/3 octave band bounds
#     fc = freqs_oto
#     flow = vcat(fc[2:end]/sqrt(2), fc[end-2]*sqrt(2))
#     fhi = vcat(fc[3]/sqrt(2), fc[1:end-1]*sqrt(2))
#
#     # Convert pseudo-density to amplitude
#     # dfs = vcat(freqs_oto[2]-freqs_oto[1], [freqs_oto[i+1] - freqs_oto[i] for i in 1:length(freqs_oto)-1])
#     # ampl_oto = sqrt.(Gxx_oto .* dfs)
#     # dfs = [freqs_oto[i+1] - freqs_oto[i] for i in 1:length(freqs_oto)-1]
#     # ampl_oto = sqrt.(Gxx_oto[1:end-1] .* dfs)
#
#     dfs = fhi .- flow
#     ampl_oto = sqrt.(Gxx_oto .* dfs)
#
#     # Convert amplitude to SPL
#     spls_oto = 20*log10.(ampl_oto./pref)
#
#     # return freqs_oto[1:end-1], spls_oto
#     return freqs_oto, spls_oto
# end


# """
# Convert spectral SPL from narrow band to one-third octave band.
# """
# function onethird2narrow_spl_DEPRE(freqs_oto, df, spls_oto; pref=20e-6)
#
#
#     # Convert SPL to OTO amplitude
#     ampl_oto = pref * 10 .^ (spls_oto ./ 20)
#
#     # 1/3 octave band bounds
#     fc = freqs_oto
#     flow = vcat(fc[2:end]/sqrt(2), fc[end-2]*sqrt(2))
#     fhi = vcat(fc[3]/sqrt(2), fc[1:end-1]*sqrt(2))
#
#     # Convert amplitude to pseudo-density
#     dfs_oto = fhi .- flow
#     Gxx_oto = ampl_oto.^2 ./ dfs_oto
#
#     # Convert one-third octave band to narrow band
#     freqs, Gxx = onethird2narrow_Gxx(freqs_oto, Gxx_oto, df)
#
#     # Convert pseudo-density to amplitude
#     ampl = sqrt.(Gxx .* df)
#
#     # Convert amplitude to SPL
#     spls = 20*log10.(ampl ./ pref)
#
#     return freqs, spls
# end
#
#
# """
# Convert the spectral density Gxx from one-third octave band to narrow band.
# """
# function onethird2narrow_Gxx(freqs_oto, Gxx_oto, df)
#
#     fc = freqs_oto
#     flow = vcat(fc[2:end]/sqrt(2), fc[end-2]*sqrt(2))
#     fhi = vcat(fc[3]/sqrt(2), fc[1:end-1]*sqrt(2))
#
#     freqs = collect(flow[1]:df:fhi[end])
#     Gxx = zeros(eltype(Gxx_oto), length(freqs))
#
#     for otoi in 1:length(fc)
#
#         # Calculate average spectral density in this interval
#         df_oto = fhi[otoi]-flow[otoi]
#         this_Gxx = Gxx_oto[otoi]/df_oto
#         # this_Gxx = Gxx_oto[otoi]
#         # this_Gxx = Gxx_oto[otoi]*df_oto/df
#         # this_Gxx = Gxx_oto[otoi]*df/df_oto
#
#         # Find corresponding narrow band interval
#         loi = findfirst(f -> f>=flow[otoi], freqs)
#         upi = findlast(f -> f<fhi[otoi], freqs)
#         upi = upi==nothing ? length(freqs) : upi
#
#         # Distribute average spectral density into the narrow band bins
#         Gxx[loi:upi] .= this_Gxx
#
#     end
#
#     return freqs, Gxx
# end


"""
Convert spectral SPL from narrow band to one-third octave band.
"""
function onethird2narrow_spl(freqs_oto, df, spls_oto; smooth=true)

    # Find interval of frequencies
    fc_min = freqs_oto[findfirst(spl -> spl>-55, spls_oto)]
    fc_max = freqs_oto[findlast(spl -> spl>-55, spls_oto)]
    flow_mini = findlast(f -> f<=fc_min, flow_oto)
    fhi_maxi = findfirst(f -> f>fc_max, fhi_oto)

    # # Chop the bins that are smaller than the narrow band frequency
    # # This leads to losing from small energy, but makes things simpler
    # while fhi_oto[flow_mini] - flow_oto[flow_mini] < df && flow_mini < length(fc_oto)
    #     flow_mini += 1
    # end

    flow_min = flow_oto[flow_mini]
    fhi_max = fhi_oto[fhi_maxi]

    freqs = collect(flow_min:df:fhi_max)
    spls = zeros(eltype(spls_oto), length(freqs))

    for otoi in 1:length(freqs_oto)

        # Find corresponding narrow band interval
        loi = findfirst(f -> f>=flow_oto[flow_mini + otoi], freqs)
        upi = findlast(f -> f<fhi_oto[flow_mini + otoi], freqs)
        loi = loi==nothing ? length(freqs) : loi
        upi = upi==nothing ? 1 : upi

        # println("$otoi\t$(flow_oto[otoi])\t$(fhi_oto[otoi])\t$loi\t$upi\t$(freqs[loi])\t$(freqs[upi])")
        # println("$(upi-loi+1)\t$(log10(upi-loi+1))\t$(spls_oto[otoi] - log10(upi-loi+1))")

        # Distribute this OTO spl bin into the narrow band SPL bins
        spls[loi:upi] .+= spls_oto[otoi] - 10*log10(upi-loi+1)

    end

    # Smooth curve up through spline if requested
    if smooth
        spls_oto_dist = zeros(eltype(spls_oto), length(spls_oto))

        # Find distributed OTO SPL values
        for (i, freq_oto) in enumerate(freqs_oto)

            # Range where this OTO freq could be
            rng = searchsorted(freqs, freq_oto)
            rng = length(rng) >= 1 ? rng : rng.stop:rng.start

            # Pick the first valid position in that range
            freq_i = rng[findfirst( i -> 1 <= i <= length(spls), rng)]

            spls_oto_dist[i] = spls[freq_i]
        end

        # Interpolate the distributed OTO SPL through spline
        spls[:] .= math.akima(freqs_oto, spls_oto_dist, freqs)
    end

    # Bring empty bins down to minimum dB
    for i in 1:length(spls)
        if spls[i] == 0 || spls[i] < -55
            spls[i] = -55
        end
    end


    return freqs, spls
end
