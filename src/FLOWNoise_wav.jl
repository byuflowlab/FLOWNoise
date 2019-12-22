#=##############################################################################
# DESCRIPTION
    Creation and processing of WAV audio files.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################

"""
# ARGUMENTS
* `in_samplerate::Real`         : Sample rate of the input pressure data WAV
                                    file in Hz.

# OPTIONAL ARGUMENTS
* `out_samplerate::Real`        : Sample rate of the output WAV file in Hz.
* `secs::Real`                  : Total length of WAV file in seconds.
* `usefreqs::Real`              : Frequencies to use (see `sampleiFFT`).
"""
function generate_WAV(pressuredata::Array{T, 1}, in_samplerate::Real,
                        filename::String; path="", out_samplerate::Real=32000,
                        secs::Real=4.0, usefreqs::Real=0.75, plot_verif=false
                            ) where {T<:Real}

    # Calculate FFT coefficients
    fftdata = fft(pressuredata)
    # freqs = in_samplerate * (0:size(fftdata, 1)-1)/size(fftdata, 1)

    samples = out_samplerate*secs      # Number of samples on WAV file

    # Generate sound wave
    y = [sampleiFFT(i/out_samplerate, fftdata, in_samplerate;
                                    usefreqs=usefreqs) for i in 0:(samples-1)]

    # Normalizes the sound wave
    maxy = max(abs(maximum(y)), abs(minimum(y)))
    y[:] .= y/maxy

    # Write wave to WAV file
    WAV.wavwrite(y, joinpath(path, filename), Fs=out_samplerate)

    if plot_verif
        period = length(pressuredata)/in_samplerate  # Period of input wave
        maxt = 4*period                             # Maximum time to plot
        maxi = ceil(Int, maxt*out_samplerate)       # Max index to plot
        plt.plot((0:length(pressuredata)-1)/in_samplerate, pressuredata, "-r",
                                                                  label="Input")
        plt.plot((0:maxi-1)/out_samplerate-period, maxy*y[1:maxi], ",k",
                                                              label="Resampled")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=false)
    end

    return y, out_samplerate
end


"""
    `sampleiFFT(t, fftdata, fftsamplerate; usefreqs=1/3)`

Given the output `fftdata = FFTW.fft(data)` and the samplerate of `data`,
it will sample the inverse Fourier transform at the time `t` and return
the value in time domain.

`usefreqs` indicates what portion of the fft frequencies to sample, with
0 indicating none, 0.5. indicating the first half of them, and 1
indicating all of them.
"""
function sampleiFFT(t::Real, fftdata::Array{C, 1}, fftsamplerate::Real;
                        usefreqs::Real=1/3) where {C<:Complex}
    val = 0

    for (i, amp) in enumerate(fftdata)

        # This condition avoids including under-calculated freqs
        if i >= length(fftdata)*usefreqs
            break
        end

        # Sample this frequency
        freq = fftsamplerate*(i-1)/length(fftdata)
        val += real(amp)*cos(2*pi*t*freq) - imag(amp)*sin(2*pi*t*freq)

    end

    # NOTE: I'm not sure why I need to multiply by 2 here, but if I don't then
    # I get half the amplitude of the input
    return 2*val/length(fftdata)
end
