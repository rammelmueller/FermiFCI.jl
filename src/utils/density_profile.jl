function compute_density_profile(::Type{T}, xgrid::Array{DType,1}, obdm::Array{DType,2})::Array{DType,1} where T<:Orbital
    """ Takes the one-body density matrix and returns the spatial density-profile
        on the specified position grid.
    """
    density = zeros(DType, (length(xgrid),))
    for i=1:size(obdm)[1]
        oi = T(i)
        di = conj.(oi.(xgrid)) 
        for j=1:size(obdm)[2]
            oj = T(j)
            density .+= di.* oj.(xgrid) .* obdm[i,j]
        end
    end
    return -density
end

#
# def get_density_profile_1D(self, x_grid, flavor, normalize=True):
#     """ Returns the density-profile for the given coordinate grid.
#     """
#     density = np.zeros_like(x_grid)
#     # for k, x in enumerate(x_grid):
#     for a in range(self.n_basis):
#         for b in range(self.n_basis):
#             density += np.conj(self.sp_orbital(a,x_grid)) * self.sp_orbital(b,x_grid) * self.obdm[flavor][a,b]
#
#     if normalize:
#         norm = np.sum(density) / self.n_part[flavor] * (x_grid[1]-x_grid[0])
#         density /= norm
#
#     return density
