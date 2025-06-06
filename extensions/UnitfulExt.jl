module UnitfulExt

using AstroCoords
using Unitful

# Remove distance from CartesianD coordinates (convert to unit vectors)
function AstroCoords.remove_distance(coord::Coordinate{<:AbstractFrame, CartesianD{T}}) where {T <: Unitful.Quantity}
    return Coordinate(coord.frame, 
        Cartesian(upreffer.(coord.representation.v))  # Cartesian constructor normalizes automatically
    )
end


end