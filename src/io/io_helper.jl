#===============================================================================

    io_helper.jl - LR, June 2021

    Some useful functionality for nicer output.

===============================================================================#
struct MemoryTag
    """ Human readable memory tag.
    """
    bytes::Integer
end

function Base.show(io::IO, tag::MemoryTag)
    print(io, round(tag.bytes/1024^2, digits=3), "MB")
end
