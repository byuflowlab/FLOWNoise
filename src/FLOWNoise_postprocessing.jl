#=##############################################################################
# DESCRIPTION
    Tools for post-processing, plotting, and visualizing acoustic outputs.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2021
  * License   : MIT
=###############################################################################

"""
    `read_data(dataset_infos; datasets_psw=def_datasets_psw,
                datasets_bpm=def_datasets_bpm,
                psw_fieldnames=["pressure", "spl_spectrum", "OASPLdB", "OASPLdBA"],
                bpm_fieldnames=["spl_spectrum", "splA_spectrum", "OASPLdB", "OASPLdBA"],
                optargs...)`

Read a series of datasets comprised of PSU-WOPWOP and BPM outputs and stores
them under the dictionaries `datasets_psw` and `datasets_bpm`.

# ARGUMENTS
* `dataset_infos`    : Array of tuples (see example).

# OPTIONAL ARGUMENTS
* `psw_fieldnames`   : Fields to read from PSU-WOPWOP.
* `bpm_fieldnames`   : Fields to read from BPM.
* `optargs...`       : Optional arguments passed to `fetch_pswdataset(...)` and
                        `fetch_bpmdataset(...)`.

# EXAMPLE

```julia
# Dataset to read and associated information
dataset_infos = [   # (label, PSW solution, BPM solution, ...)

                    ("FLOWUnsteady w/BEMT",
                        "temps/example_dji9443_00_psw00/runcase/",
                        "temps/example_dji9443_00_bpm00",
                        "-", "r"),

                    ("FLOWUnsteady w/VPM",
                        "temps/example_dji9443_01_psw00/runcase/",
                        "temps/example_dji9443_01_bpm00",
                        "-", "b")
                ]

datasets_psw = Dict()     # Stores PSW data in this dictionary
datasets_bpm = Dict()     # Stores BPM data in this dictionary

# Read datasets and stores them in dictionaries
read_data(dataset_infos; datasets_psw=datasets_psw, datasets_bpm=datasets_bpm)
```
"""
function read_data(dataset_infos;
                    datasets_psw=def_datasets_psw, datasets_bpm=def_datasets_bpm,
                    psw_fieldnames=["pressure", "spl_spectrum", "spl_octFilt_spectrum", "OASPLdB", "OASPLdBA"],
                    bpm_fieldnames=["spl_spectrum", "splA_spectrum", "OASPLdB", "OASPLdBA", "frequencies"],
                    optargs...)

    # Read each dataset
    for (lbl, psw_read_path, bpm_read_path) in dataset_infos

        fetch_pswdataset(psw_read_path; datasets=datasets_psw, fieldnames=psw_fieldnames, optargs...)
        fetch_bpmdataset(bpm_read_path; datasets=datasets_bpm, fieldnames=bpm_fieldnames, optargs...)

    end
end



"""
`plot_pressure(dataset_infos, microphones, RPM, sph_ntht, pangle;
fieldname="pressure", datasets_psw=def_datasets_psw)`

Plots the pressure waveform of the tonal component of noise calculated by
PSU-WOPWOP at the requested microphone positions of a 2D grid.

# ARGUMENTS
* `dataset_infos`     : Array of tuples containing datasets to plot (see
                            example).
* `microphones`       : Array of microphone positions to plot.
* `RPM`               : RPM of the simulation (used to scale the time).
* `sph_ntht`          : Number of azimuthal nodes in the grid of observers (used
                            to convert `microphones` into node indices).
* `pangle::Function`  : Function to convert from indices of grid nodes to
                            degrees (or any other metric represented by the
                            azimuthal dimension of the grid).

# OPTIONAL ARGUMENTS
* `datasets_psw`      : Datasets where to read data from.

# EXAMPLE

```julia
# Dataset to read and associated information
dataset_infos = [   # (label, PSW solution, BPM solution, line style, color)

                    ("FLOWUnsteady w/BEMT",
                        "temps/example_dji9443_00_psw00/runcase/",
                        "temps/example_dji9443_00_bpm00",
                        "-", "r"),

                    ("FLOWUnsteady w/VPM",
                        "temps/example_dji9443_01_psw00/runcase/",
                        "temps/example_dji9443_01_bpm00",
                        "-", "b")
                ]

datasets_psw = Dict()     # Stores PSW data in this dictionary
datasets_bpm = Dict()     # Stores BPM data in this dictionary

# Read datasets and stores them in dictionaries
read_data(dataset_infos; datasets_psw=datasets_psw, datasets_bpm=datasets_bpm)

# Plot pressure waveform
microphones  = [-45, -90]            # (deg) microphones to plot
RPM          = 5400                  # RPM of solution

plot_pressure(dataset_infos, microphones, RPM, sph_ntht, pangle; datasets_psw=datasets_psw)
```
"""
function plot_pressure(dataset_infos, microphones, RPM, sph_ntht, pangle::Function;
                                fieldname="pressure", datasets_psw=def_datasets_psw,
                                xlims=[0, 3])
    # microphones  = [-45, -90]            # (deg) microphones to plot
    # RPM          = 5400                  # RPM of solution

    # fieldname = "pressure"               # Field to plot
    hash = Dict((round(pangle(mici), digits=1), mici) for mici in 1:sph_ntht)
    # mics         = Int.((-microphones .+ 180) * sph_ntht/360 .+ 1)   # Index of every microphone
    mics         = [hash[deg] for deg in microphones]

    for mici in mics

        plt.figure("$(fieldname)-$(mici)", figsize=[7*1*1.4, 3*3]*0.7)
        plt.suptitle("$(Int(round(pangle(mici))))"*L"$^\circ$ microphone")

        for (ploti, (ylbl, dlbl)) in enumerate([
                                                ("Total pressure (Pa)", "TotalAcousticPressure"),
                                                ("Loading pressure (Pa)", "LoadingAcousticPressure"),
                                                ("Thickness pressure (Pa)", "ThicknessAcousticPressure"),
                                              ])
            plt.subplot(310+ploti)
            if ploti==3; plt.xlabel("Rotor revolution"); end;
            plt.xlim(xlims)
            plt.ylabel(ylbl)
            plt.grid(true, which="major", color="0.8", linestyle="--")
            plt.grid(true, which="minor", color="0.8", linestyle="--")


            for (di, (lbl, read_path, _, stl, clr)) in enumerate(dataset_infos)

                data = datasets_psw[read_path][fieldname]
                xi = data["hs"]["ObserverTimes"]
                yi = data["hs"][dlbl]

                plt.plot(data["field"][mici, 1, 2:end, xi]*(RPM/60),
                         data["field"][mici, 1, 2:end, yi],  stl, alpha=0.8, label=lbl, color=clr)

            end

            if ploti==1; plt.legend(loc="center left", bbox_to_anchor=(1, 0.5)); end;
            # if ploti==1; plt.legend(loc="best"); end;
        end

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    end
end

"""
    `plot_spectrum_spl(dataset_infos, microphones, BPF, sph_ntht, pangle::Function;
                                datasets_psw=def_datasets_psw,
                                datasets_bpm=def_datasets_bpm,
                                onethirdoctave=false,
                                Aweighted=false,
                                plot_csv=nothing, csv_lbl="Experimental", csv_stl="k", csv_weight=true,
                                xBPF=true, BPF_lines=0,
                                xlims=[6*10.0^(-1.0), 1.75e2], ylims=[0, 60],
                                Pref=L" re \$20\\mu\$Pa"
                                )`

Plot the SPL spectrum in `dataset_infos` at the given microphones.

# ARGUMENTS
* `dataset_infos`     : Array of tuples containing datasets to plot (see
                            example).
* `microphones`       : Array of microphone positions to plot.
* `BPF`               : Blade passing frequency (used to scale the frequency).
* `sph_ntht`          : Number of azimuthal nodes in the grid of observers (used
                            to convert `microphones` into node indices).
* `pangle::Function`  : Function to convert from indices of grid nodes to
                            degrees (or any other metric represented by the
                            azimuthal dimension of the grid).

# OPTIONAL ARGUMENTS
* `datasets_psw`      : PSU-WOPWOP datasets where to read data from.
* `datasets_bpm`      : BPM datasets where to read data from.
* `onethirdoctave`    : Converts the plot to 1/3 octave bands, otherwise it'll
                            convert BPM's 1/3 bands to narrow bands of PSW.
* `Aweighted`         : Plot A-weighted SPL.
* `plot_csv`          : Give it the path to a CSV file to include in the plot.
* `xBPF`              : Scale frequencies by BPF if true.
* `BPF_lines`         : Number of lines of BPF multiples to draw in the plot.

# EXAMPLE
```julia
# Dataset to read and associated information
dataset_infos = [   # (label, PSW solution, BPM solution, BPM freq bins, line style, color)
                    ("FLOWUnsteady w/BEMT",
                        "temps/example_dji9443_00_psw00/runcase/",
                        "temps/example_dji9443_00_bpm00",
                        "-", "r")
                ]

datasets_psw = Dict()     # Stores PSW data in this dictionary
datasets_bpm = Dict()     # Stores BPM data in this dictionary

# Read datasets and stores them in dictionaries
FLOWNoise.read_data(dataset_infos; datasets_psw=datasets_psw, datasets_bpm=datasets_bpm)

microphones  = [-45]                 # (deg) microphone to plot
Aweighted    = false                 # Plot A-weighted SPL

# Experimental data from Zawodny et al., Fig. 9
filename = joinpath("data", "zawodnydigitization2", "zawodny_dji9443_spl_5400_01.csv")
plot_experimental = [(filename, "Experimental", "k", Aweighted, [])]

# Plot spectrum in Hz with blade-passing-frequency lines
FLOWNoise.plot_spectrum_spl(dataset_infos, microphones, BPF, sph_ntht, pangle;
                              datasets_psw=datasets_psw, datasets_bpm=datasets_bpm,
                              Aweighted=Aweighted,
                              plot_csv=plot_experimental,
                              xBPF=false, xlims=[100, 3e4], BPF_lines=21)
```

"""
function plot_spectrum_spl(dataset_infos, microphones, BPF, sph_ntht, pangle::Function;
                            datasets_psw=def_datasets_psw,
                            datasets_bpm=def_datasets_bpm,
                            onethirdoctave=false,
                            psworg_onethirdoctave=false,
                            smooth_narrow=true,
                            Aweighted=false,
                            components=false,
                            add_broadband=true,
                            plot_csv=[],
                            xBPF=true, BPF_lines=0,
                            figname="automatic",
                            xlims=[6*10.0^(-1.0), 1.75e2], ylims=[0, 60],
                            Pref=L" re $20\mu$Pa",
                            verbose=true
                            )

    # microphone   = -45                   # (deg) microphone to plot
    # RPM          = 5400                  # RPM of solution
    # BPF          = 2*RPM/60              # Blade passing frequency

    fieldname_psw = onethirdoctave && psworg_onethirdoctave ? "spl_octFilt_spectrum" : "spl_spectrum"    # Field to plot
    fieldname_bpm = "spl"*"A"^Aweighted*"_spectrum"        # Field to plot


    hash = Dict((round(pangle(mici), digits=1), mici) for mici in 1:sph_ntht)
    # micis         = Int.((-microphones .+ 180) * sph_ntht/360 .+ 1)    # Index of each microphone
    micis         = [hash[deg] for deg in microphones]

    xscaling = (1/BPF)^xBPF

    aux = "SPL"*(onethirdoctave ? L"_{1/3}" : "")
    ylbl = Aweighted ? "A-weighted "*aux*"\n(dBA"*Pref*")" : aux*" (dB"*Pref*")"

    for mici in micis

        _figname = figname=="automatic" ? "$(fieldname_bpm)-$(mici)" : figname
        plt.figure(_figname, figsize=[7, 3.5])

        for (csv_filename, lbl, stl, weight, optargs) in plot_csv
            data = CSV.read(csv_filename, DataFrames.DataFrame, datarow=1)
            spl_exp = weight ? aWeight.(data[!, 1], data[!, 2]) : data[!, 2]
            plt.plot(data[!, 1]*xscaling, spl_exp, stl; label=lbl, optargs...)
        end

        # Plot datasets
        for (lbl, read_path_psw, read_path_bpm, stl, clr) in dataset_infos

            # Fetch tonal noise
            data_psw = datasets_psw[read_path_psw][fieldname_psw]
            xi = data_psw["hs"]["Frequency"]
            yi = data_psw["hs"]["Total_dB"*"A"^Aweighted]
            freqs_psw = data_psw["field"][mici, 1, 2:end, xi]
            org_spl_psw = data_psw["field"][mici, 1, 2:end, yi]

            # Fetch broadband noise
            data_bpm = datasets_bpm[read_path_bpm][fieldname_bpm]
            freqs_bpm = datasets_bpm[read_path_bpm]["frequencies"]["field"][:, 1]
            org_spl_bpm = data_bpm["field"][:, mici]

            # Convert narrow -> 1/3 or 1/3 -> narrow band to make tonal and
            # broadband comparable
            if onethirdoctave
                if !psworg_onethirdoctave
                    freqs_psw, org_spl_psw = narrow2onethird_spl(freqs_psw, org_spl_psw)
                end
            else
                df = (freqs_psw[2]-freqs_psw[1])
                println("Converting BPM to narrow band of bin width $(round(df, digits=1)) Hz")
                freqs_bpm, org_spl_bpm = onethird2narrow_spl(freqs_bpm, df, org_spl_bpm; smooth=smooth_narrow)
            end

            # Define the range of frequency as the union of both components
            freqs = freqs_psw                       # Grab PSW frequency range
            if freqs_psw[end] < freqs_bpm[end]  # Add BPM range
                bpm_i = findfirst( f -> f > freqs_psw[end], freqs_bpm)
                freqs = vcat(freqs, freqs_bpm[bpm_i:end])
                spl_psw = vcat(org_spl_psw, [-30.0 for i in bpm_i:length(freqs_bpm)])
            else
                spl_psw = org_spl_psw
            end

            # Interpolate broadband data into the same frequencies than tonal data
            spl_bpm = math.akima(freqs_bpm, org_spl_bpm, freqs)

            # Add tonal and broadband SPL together
            spl = addSPL(spl_psw, spl_bpm)

            if components
                plt.plot(freqs*xscaling, spl, "-", alpha=0.75, label=lbl*" --- Total", color=clr)
                plt.plot(freqs_bpm*xscaling, org_spl_bpm, ":", alpha=0.75, label=lbl*" --- BPM", color=clr)
                plt.plot(freqs_psw*xscaling, org_spl_psw, "--", alpha=0.75, label=lbl*" --- FW-H", color=clr)
            elseif !add_broadband
                plt.plot(freqs_psw*xscaling, org_spl_psw, stl, alpha=0.75, label=lbl, color=clr)
            else
                plt.plot(freqs*xscaling, spl, ("o"^onethirdoctave)*stl, alpha=0.75, label=lbl, color=clr)
            end

        end

        if BPF_lines != 0
            for BPFi in 1:BPF_lines
                this_x = BPFi*BPF*xscaling
                plt.plot(this_x*ones(100), range(ylims[1], ylims[2], length=100),
                                                               ",k", alpha=0.25)
                if BPFi<=8
                    plt.annotate("$BPFi"*(BPFi==1 ? L"$^\mathrm{st}$" :
                                          BPFi==2 ? L"$^\mathrm{nd}$" :
                                          BPFi==3 ? L"$^\mathrm{rd}$" : L"$^\mathrm{th}$"),
                                  (this_x, ylims[2]), (this_x, (1-0.1-0.001*BPFi^2.35)*(ylims[2]-ylims[1])), alpha=0.25)
                end
            end
        end

        plt.title("$(Int(round(pangle(mici))))"*L"$^\circ$ microphone")
        plt.xlim(xlims)
        plt.xscale("log")
        plt.xlabel("Frequency"*(xBPF ? " / BPF" : " (Hz)"))
        plt.ylim(ylims)
        plt.ylabel(ylbl)
        if BPF_lines==0 && BPF_lines!=-1
            plt.grid(true, which="major", color="0.8", linestyle=":")
            plt.grid(true, which="minor", color="0.8", linestyle=":")
        end
        plt.legend(loc="best", frameon=false)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    end
end

"""
    `plot_directivity_splbpf(dataset_infos, BPFi, BPF, pangle::Function;
                                        datasets_psw=def_datasets_psw,
                                        datasets_bpm=def_datasets_bpm,
                                        plot_csv=[],
                                        thetalims=180*[1, -1],
                                        thetagrids=collect(-180:45/2:180)[[3, 5, 7, 8, 9, 10, 11, 13, 15, 17]],
                                        rlbl_pos=90,
                                        rticks=40:4:52, rlims=[40, 54], rorigin=36,
                                        Pref=L" re \$20\\mu\$Pa")`

Plot the directivity of the requested multiple `BPFi` of the blade-passing
frequency `BPF`.

# ARGUMENTS
* `dataset_infos`     : Array of tuples containing datasets to plot (see
                            example).
* `pangle::Function`  : Function to convert from indices of grid nodes to
                            degrees (or any other metric represented by the
                            azimuthal dimension of the grid) as `pangle(i)`.

# OPTIONAL ARGUMENTS
* `datasets_psw`      : PSU-WOPWOP datasets where to read data from.
* `datasets_bpm`      : BPM datasets where to read data from.
* `plot_csv`          : Give it the path to a CSV file to include in the plot.

# EXAMPLE
```julia

# Dataset to read and associated information
dataset_infos = [   # (label, PSW solution, BPM solution, BPM freq bins, line style, color)
                    ("FLOWUnsteady w/BEMT",
                        "temps/example_dji9443_00_psw00/runcase/",
                        "temps/example_dji9443_00_bpm00",
                        "-", "r")
                ]

datasets_psw = Dict()     # Stores PSW data in this dictionary
datasets_bpm = Dict()     # Stores BPM data in this dictionary

# Read datasets and stores them in dictionaries
FLOWNoise.read_data(dataset_infos; datasets_psw=datasets_psw, datasets_bpm=datasets_bpm)

BPFi = 1                 # Multiple of blade-passing frequency to plot

# Experimental and computational data reported by Zawodny et al., Fig. 14
filename1 = joinpath("data", "zawodnydigitization2", "zawodny_fig14_topright_exp00.csv")
filename2 = joinpath("data", "zawodnydigitization2", "zawodny_fig14_topright_of00.csv")
filename3 = joinpath("data", "zawodnydigitization2", "zawodny_fig14_topright_pas00.csv")

plot_experimental = [(filename1, "Experimental", "ok", Aweighted, []),
                     (filename2, "OVERFLOW2", "-c", Aweighted, [(:alpha, 0.8)]),
                     (filename3, "BEM+ANOPP", ":b", Aweighted, [(:alpha, 0.8)])]

# Plot SPL directivity of first blade-passing frequency
FLOWNoise.plot_directivity_splbpf(dataset_infos, BPFi, BPF, pangle;
                                    datasets_psw=datasets_psw,
                                    datasets_bpm=datasets_bpm,
                                    plot_csv=plot_experimental,
                                    rticks=40:4:52, rlims=[40, 54], rorigin=36)
```
"""
function plot_directivity_splbpf(dataset_infos, BPFi, BPF, pangle::Function;
                                    datasets_psw=def_datasets_psw,
                                    datasets_bpm=def_datasets_bpm,
                                    plot_csv=[],
                                    thetalims=180*[1, -1],
                                    thetagrids=collect(-180:45/2:180)[[3, 5, 7, 8, 9, 10, 11, 13, 15, 17]],
                                    rlbl_pos=90,
                                    rticks=40:4:52, rlims=[40, 54], rorigin=36,
                                    Pref=L" re $20\mu$Pa")

    # BPFi         = 1                     # BPF multiple to plot
    fieldname_psw = "spl_spectrum"         # Field to plot
    fieldname_bpm = "spl_spectrum"         # Field to plot

    plt.figure("$(fieldname_bpm)-$(BPFi)", figsize=[3, 3]*2)

    for (csv_filename, lbl, stl, weight, optargs) in plot_csv
        data = CSV.read(csv_filename, DataFrames.DataFrame, datarow=1)
        plt.polar(pi/180*data[!, 2], data[!, 1], stl; label=lbl, optargs...)
    end

    # Plot datasets
    for (lbl, read_path_psw, read_path_bpm, stl, clr) in dataset_infos

        # Fetch tonal noise
        data_psw = datasets_psw[read_path_psw][fieldname_psw]
        yi = data_psw["hs"]["Total_dB"]
        fi = data_psw["hs"]["Frequency"]

        df = data_psw["field"][1, 1, 2, fi] - data_psw["field"][1, 1, 1, fi] # Frequency step
        freqi = ceil(Int, BPFi*BPF/df + 1)                           # Frequency index
        freq = data_psw["field"][1, 1, freqi, fi]                    # Frequency
        # _lbl = " @ $(Int(round(freq, digits=0))) Hz"                 # Frequency string
        _lbl = ""

        spl_psw = data_psw["field"][:, 1, freqi, yi]

        # Fetch broadband noise
        data_bpm = datasets_bpm[read_path_bpm][fieldname_bpm]
        freqs_bpm = datasets_bpm[read_path_bpm]["frequencies"]["field"][:, 1]
        spl_bpm = [math.akima(freqs_bpm, data_bpm["field"][:, mici], freq) for mici in 1:length(spl_psw)]

        # Add tonal and broadband SPL together
        # spl = addSPL(spl_psw, spl_bpm)
        # NOTE: Here we ignore the broadband component since it is in a 1/3 octave band
        spl = spl_psw

        angles = pangle.(1:length(spl))
        plt.polar(pi/180*angles, spl, stl, label=lbl*_lbl, alpha=0.8, color=clr)
    end

    ax = plt.gca()
    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)
    ax.set_thetalim(pi/180*thetalims)
    ax.set_rlabel_position(rlbl_pos)
    ax.set_rticks(rticks)
    ax.set_rlim(rlims)
    ax.set_rorigin(rorigin)
    ax.set_thetagrids(thetagrids)
    ax.grid(true)
    plt.xlabel("SPL (dB"*Pref*")\n @ $(BPFi)"*(BPFi==1 ? L"$^\mathrm{st}$" :
                           BPFi==2 ? L"$^\mathrm{nd}$" :
                           BPFi==3 ? L"$^\mathrm{rd}$" : L"$^\mathrm{th}$")*" BPF ($(ceil(Int, BPFi*BPF)) hz)")
    plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))
    plt.tight_layout();
end

"""
    `plot_directivity_oaspl(dataset_infos, pangle::Function;
                                        datasets_psw=def_datasets_psw,
                                        datasets_bpm=def_datasets_bpm,
                                        Aweighted=false,
                                        plot_csv=[],
                                        thetalims=180*[1, -1],
                                        thetagrids=collect(-180:45/2:180)[[3, 5, 7, 8, 9, 10, 11, 13, 15, 17]],
                                        rlbl_pos=90,
                                        rticks=40:10:70, rlims=[40, 72], rorigin=30,
                                        Pref=L" re \$20\\mu\$Pa")`

Plot overall SPL directivity.

# ARGUMENTS
* `dataset_infos`     : Array of tuples containing datasets to plot (see
                            example).
* `pangle::Function`  : Function to convert from indices of grid nodes to
                            degrees (or any other metric represented by the
                            azimuthal dimension of the grid) as `pangle(i)`.

# OPTIONAL ARGUMENTS
* `datasets_psw`      : PSU-WOPWOP datasets where to read data from.
* `datasets_bpm`      : BPM datasets where to read data from.
* `Aweighted`         : Plot A-weighted SPL.
* `plot_csv`          : Give it the path to a CSV file to include in the plot.

# EXAMPLE
```julia

# Dataset to read and associated information
dataset_infos = [   # (label, PSW solution, BPM solution, BPM freq bins, line style, color)
                    ("FLOWUnsteady w/BEMT",
                        "temps/example_dji9443_00_psw00/runcase/",
                        "temps/example_dji9443_00_bpm00",
                        "-", "r")
                ]

datasets_psw = Dict()     # Stores PSW data in this dictionary
datasets_bpm = Dict()     # Stores BPM data in this dictionary

# Read datasets and stores them in dictionaries
FLOWNoise.read_data(dataset_infos; datasets_psw=datasets_psw, datasets_bpm=datasets_bpm)

Aweighted = false

# Experimental and computational data reported by Zawodny et al., Fig. 12
filename1 = joinpath("data", "zawodnydigitization2", "zawodny_fig12_left_5400_00.csv")

plot_experimental = [(filename1, "Experimental", "ok", Aweighted, [])]

# Plot OASPL directivity
FLOWNoise.plot_directivity_oaspl(dataset_infos, pangle;
                                    datasets_psw=datasets_psw,
                                    datasets_bpm=datasets_bpm,
                                    Aweighted=Aweighted,
                                    plot_csv=plot_experimental,
                                    rticks=40:10:70, rlims=[40, 72], rorigin=30)
```
"""
function plot_directivity_oaspl(dataset_infos, pangle::Function;
                                    datasets_psw=def_datasets_psw,
                                    datasets_bpm=def_datasets_bpm,
                                    Aweighted=false,
                                    plot_csv=[],
                                    thetalims=180*[1, -1],
                                    thetagrids=collect(-180:45/2:180)[[3, 5, 7, 8, 9, 10, 11, 13, 15, 17]],
                                    rlbl_pos=90,
                                    rticks=40:10:70, rlims=[40, 72], rorigin=30,
                                    Pref=L" re $20\mu$Pa")

    fieldname_psw = "OASPLdB"*"A"^Aweighted # Field to plot
    fieldname_bpm = fieldname_psw           # Field to plot

    plt.figure("$(fieldname_bpm)", figsize=[3, 3]*2)

    for (csv_filename, lbl, stl, weight, optargs) in plot_csv
        data = CSV.read(csv_filename, DataFrames.DataFrame, datarow=1)
        plt.polar(pi/180*data[!, 2], data[!, 1], stl; label=lbl, optargs...)
    end

    # Plot datasets
    for (lbl, read_path_psw, read_path_bpm, stl, clr) in dataset_infos

        # Fetch tonal noise
        data_psw = datasets_psw[read_path_psw][fieldname_psw]
        yi = data_psw["hs"]["TotalOASPLdB"*"A"^Aweighted]
        oaspl_psw = data_psw["field"][:, 1, 1, yi]

        # Fetch broadband noise
        data_bpm = datasets_bpm[read_path_bpm][fieldname_bpm]
        oaspl_bpm = data_bpm["field"][:]

        # Add tonal and broadband OASPL together
        oaspl = addSPL(oaspl_psw, oaspl_bpm)

        angles = pangle.(1:length(oaspl))
        plt.polar(pi/180*angles, oaspl, stl, label=lbl, alpha=0.8, color=clr)
    end

    ax = plt.gca()
    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)
    ax.set_thetalim(pi/180*thetalims)
    ax.set_rlabel_position(rlbl_pos)
    ax.set_rticks(rticks)
    ax.set_rlim(rlims)
    ax.set_rorigin(rorigin)
    ax.set_thetagrids(thetagrids)
    ax.grid(true)
    plt.xlabel("A-weighted "^Aweighted*"OASPL\n"*"(dB"*"A"^Aweighted*"$(Pref))")
    plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))
    plt.tight_layout();
end
