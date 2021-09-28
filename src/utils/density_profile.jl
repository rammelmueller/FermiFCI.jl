function compute_density_profile(orbital::T, xgrid::Array{DType,1}, obdm::Array{DType,2}; kwargs...)::Array{DType,1} where T<:Orbital
    """ Takes the one-body density matrix and returns the spatial density-profile
        on the specified position grid (unnormalized).
    """
    density = zeros(DType, (length(xgrid),))
    for i=1:size(obdm)[1]
        di = conj.(orbital.(i, xgrid))
        for j=1:size(obdm)[2]
            density .+= di.* orbital.(j, xgrid) .* obdm[i,j]
        end
    end
    return -density
end
