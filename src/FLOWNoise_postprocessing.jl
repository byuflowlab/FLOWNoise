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
                    psw_fieldnames=["pressure", "spl_spectrum", "OASPLdB", "OASPLdBA"],
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
                                fieldname="pressure", datasets_psw=def_datasets_psw)
    # microphones  = [-45, -90]            # (deg) microphones to plot
    # RPM          = 5400                  # RPM of solution

    # fieldname = "pressure"               # Field to plot
    mics         = Int.((-microphones .+ 180) * sph_ntht/360 .+ 1)   # Index of every microphone

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
            plt.xlim([0, 3])
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
* `Aweighted`         : Plot A-weighted SPL.
* `plot_csv`          : Give it the path to a CSV file to include in the plot.
* `xBPF`              : Scale frequencies by BPF if true.
* `BPF_lines`         : Number of lines of BPF multiples to draw in the plot.

# EXAMPLE
#
"""
function plot_spectrum_spl(dataset_infos, microphones, BPF, sph_ntht, pangle::Function;
                            datasets_psw=def_datasets_psw,
                            datasets_bpm=def_datasets_bpm,
                            Aweighted=false,
                            plot_csv=[],
                            xBPF=true, BPF_lines=0,
                            xlims=[6*10.0^(-1.0), 1.75e2], ylims=[0, 60],
                            Pref=L" re $20\mu$Pa"
                            )

    # microphone   = -45                   # (deg) microphone to plot
    # RPM          = 5400                  # RPM of solution
    # BPF          = 2*RPM/60              # Blade passing frequency

    fieldname_psw = "spl_spectrum"        # Field to plot
    fieldname_bpm = "spl"*"A"^Aweighted*"_spectrum"        # Field to plot

    micis         = Int.((-microphones .+ 180) * sph_ntht/360 .+ 1)    # Index of each microphone
    xscaling = (1/BPF)^xBPF

    ylbl = Aweighted ? "A-weighted SPL\n(dBA"*Pref*")" : "SPL (dB"*Pref*")"

    for mici in micis

        plt.figure("$(fieldname_bpm)-$(mici)", figsize=[7, 3.5])

        for (csv_filename, lbl, stl, weight) in plot_csv
            data = CSV.read(csv_filename, DataFrames.DataFrame)
            spl_exp = weight ? aWeight.(data[!, 1], data[!, 2]) : data[!, 2]
            plt.plot(data[!, 1]*xscaling, spl_exp, stl, label=lbl)
        end

        # Plot datasets
        for (lbl, read_path_psw, read_path_bpm, stl, clr) in dataset_infos

            # Fetch tonal noise
            data_psw = datasets_psw[read_path_psw][fieldname_psw]
            xi = data_psw["hs"]["Frequency"]
            yi = data_psw["hs"]["Total_dB"*"A"^Aweighted]
            freqs_psw = data_psw["field"][mici, 1, 2:end, xi]
            spl_psw = data_psw["field"][mici, 1, 2:end, yi]

            # Fetch broadband noise
            data_bpm = datasets_bpm[read_path_bpm][fieldname_bpm]
            freqs_bpm = datasets_bpm[read_path_bpm]["frequencies"]["field"][:, 1]
            org_spl_bpm = data_bpm["field"][:, mici]

            # Define the range of frequency as the union of both components
            freqs = freqs_psw                       # Grab PSW frequency range
            if freqs_psw[end] < freqs_bpm[end]  # Add BPM range
                bpm_i = findfirst( f -> f > freqs_psw[end], freqs_bpm)
                freqs = vcat(freqs, freqs_bpm[bpm_i:end])
                spl_psw = vcat(spl_psw, [-30.0 for i in bpm_i:length(freqs_bpm)])
            end

            # Interpolate broadband data into the same frequencies than tonal data
            spl_bpm = math.akima(freqs_bpm, org_spl_bpm, freqs)

            # Add tonal and broadband SPL together
            spl = addSPL(spl_psw, spl_bpm)

            plt.plot(freqs*xscaling, spl, stl, alpha=0.75, label=lbl, color=clr)

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
        if BPF_lines==0
            plt.grid(true, which="major", color="0.8", linestyle=":")
            plt.grid(true, which="minor", color="0.8", linestyle=":")
        end
        plt.legend(loc="best", frameon=false)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    end
end

# function plot_directivity_splbpf()
