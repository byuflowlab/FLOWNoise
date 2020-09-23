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
function autospec(x, fs::Real, ns::Int, overlap::Real; window=hanning, density=false)

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

    OASPL = 20*log10(sqrt(sum(Gxx))/2e-5)

    # Normalize by the frequency bin to get density
    if density
        Gxx /= delta_freq
    end

    return Gxx, freqs, OASPL
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
