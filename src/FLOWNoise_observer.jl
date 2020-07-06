#=##############################################################################
# DESCRIPTION
    Creating and postprocessing of observer positions.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez and Tyler Critchfield
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2019
  * License   : MIT
=###############################################################################

"""
    `observer_sphere(R, nR, ntht, nphi; Rmin=0.1, thtmin=5, thtmax=175)`

Returns a spherical grid of observers.

# ARGUMENTS
* `R::Real`                     : Sphere radius.
* `nR::Int`                     : Number of radial sections (give 0 for a hollow
                                    sphere).
* `ntht::Int`                   : Number of the elevation sections.
* `nphi::Int`                   : Number of azimuth sections.

# OPTIONAL ARGUMENTS
* `Rmin::Real`                  : Lower bound of R to mesh (given as a fraction
                                    of R).
* `thtmin::Real`                : Lower bound of elevation to mesh (in degrees).
* `thtmax::Real`                : Upper bound of elevation to mesh (in degrees).
* `phimin::Real`                : Lower bound of azimuth to mesh (in degrees).
* `phimax::Real`                : Upper bound of azimuth to mesh (in degrees).
* `half::Bool`                  : Half or full sphere.
* `C::Array{Real, 1}`           : Center [x,y,z] of the sphere.
* `rotation::Array{Real, 1}`    : Rotation of sphere (in degrees) about x, y,
                                    and z axes.
* `save_path::String`           : Saves the sphere to a this path.
* `file_name::String`           : Name of file when saving sphere.
* `fmt::String`                 : Format of the file when saving. Options:
                                    "plot3d", "vtk".

# NOTES
* Elevation angles go from 0 to 180, meanwhile azimuth angles go from 0 to 360.
* To generate a half sphere use `phimax=180`.
* Give `nphi=2*ntht` to get cells of aspect ratio 1.
* A circular line (instead of a spherical surface/volume) can be generated with
    `nphi=0` and `nR=0`. This will return a circular line laying on the
    zx-plane. Use `C` and `rotation` to translate a rotate the line.
"""
function observer_sphere(R::Real, nR::Int, ntht::Int, nphi::Int;
                            Rmin::Real=0.1,
                            thtmin=5, thtmax=175,
                            phimin=0, phimax=360,
                            C::Array{T1, 1}=zeros(3),
                            rotation::Array{T2, 1}=zeros(3),
                            save_path::Union{Nothing, String}=nothing,
                            file_name::String="sphere",
                            fmt::String="vtk"
                         ) where {T1<:Real, T2<:Real}



    P_min = [pi/180*thtmin, pi/180*phimin, R*Rmin] # Lower bound (theta, phi, r)
    P_max = [pi/180*thtmax, pi/180*phimax, R]   # Upper bound (theta, phi, r)
    NDIVS = [ntht, nphi, nR]                    # Number of divisions (cells)
                                                #  of (theta, phi, r)

    # Case that grid is a circular line
    if nphi==0 && nR==0
        P_min[3] = P_max[3]
        P_min[2] = P_max[2]
        loop_dim = Int(P_min[1]==0 && P_max[1]==2*pi) # Coordinate to loop (0==no looping)

    else

        if NDIVS[3]==0; P_min[3]=P_max[3]; end;
        loop_dim = 2*Int(P_min[2]==0 && P_max[2]==2*pi) # Coordinate to loop (0==no looping)

    end

    # Generates parametric (theta, phi, r) grid
    grid = gt.Grid(P_min, P_max, NDIVS, loop_dim)

    # Transforms the grid into a spherical cartesian space
    my_transform(X) = gt.spherical3D(vcat(X[3], X[1:2]))

    # Rotate and translate the sphere
    rotM = gt.rotation_matrix2(rotation...)
    my_transform2(X) = rotM*my_transform(X) + C

    # Apply transformation
    gt.transform!(grid, my_transform2)

    if save_path!=nothing
        if fmt=="plot3d"
            # Create .xyz file
            gt.generatePLOT3D(grid, file_name; path=save_path, trickwopwop=true)
        elseif fmt=="vtk"
            # Create .vtk file
            gt.save(grid, file_name; path=save_path)
        else
            error("Invalid file format $fmt.")
        end
    end

    return grid
end
