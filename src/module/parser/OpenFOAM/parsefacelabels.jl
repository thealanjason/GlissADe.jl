#=
# Parser to parse facelabels in OpenFOAM meshfiles.
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#
"""
    _parsefacelabels_openfoam(location_faceLabels, correction_factor)
Parse indices of the faces forming the terrain in a text file with OpenFoam format and stores it in `INT_TYPE`

## Arguments
- `location_faceLabels::String` - location of the facelabels file.
- `correction_factor::INT_TYPE` - remove incorrect entries
"""
function _parsefacelabels_openfoam(
    location_faceLabels::String,
    correction_factor::INT_TYPE[],
)
    s_facelabels_file = open(location_faceLabels, "r") do f_facelabels
        readlines(f_facelabels)
    end
    n_facelabels = parse(INT_TYPE[], s_facelabels_file[20])
    facelabels = Vector{INT_TYPE[]}(undef, n_facelabels)
    s_facelabels = @view s_facelabels_file[22:(22+n_facelabels-1)]
    STATS[] && println("Reading face labels...")
    @inbounds @maybe_threads Threads.nthreads == 1 || !THREADS[] for i = 1:n_facelabels
        facelabels[i] = parse(INT_TYPE[], s_facelabels[i]) + 1 + correction_factor
    end
    return facelabels
end
