include("coordinate_accessors.jl")

# Transformations between SphericalD (with distance) and CartesianD (with distance)
function spherical2cartesian(coord::Coordinate{Frame, SphericalD{T, D}}) where {Frame, T, D}
    return Coordinate(coord.frame, 
        CartesianD(
            SVector(
                dist(coord) * cos(lat(coord)) * cos(lon(coord)),
                dist(coord) * cos(lat(coord)) * sin(lon(coord)),
                dist(coord) * sin(lat(coord))
            )
        )
    )
end

# Transformations between Spherical (unit sphere) and Cartesian (unit vectors)
function spherical2cartesian(coord::Coordinate{Frame, Spherical{T}}) where {Frame, T}
    return Coordinate(coord.frame, 
        Cartesian(
            SVector(
                cos(lat(coord)) * cos(lon(coord)),
                cos(lat(coord)) * sin(lon(coord)),
                sin(lat(coord))
            )
        )
    )
end

# CartesianD (with distance) to SphericalD (with distance)
function cartesian2spherical(coord::Coordinate{Frame, CartesianD{T}}) where {Frame, T}
    return Coordinate(coord.frame, 
        SphericalD(
            atan(z_coord(coord), sqrt(x_coord(coord)^2 + y_coord(coord)^2)),
            atan(y_coord(coord), x_coord(coord)),
            sqrt(x_coord(coord)^2 + y_coord(coord)^2 + z_coord(coord)^2)
        )
    )
end

# Cartesian (unit vectors) to Spherical (unit sphere)
function cartesian2spherical(coord::Coordinate{Frame, Cartesian{T}}) where {Frame, T}
    return Coordinate(coord.frame, 
        Spherical(
            atan(z_coord(coord), sqrt(x_coord(coord)^2 + y_coord(coord)^2)),
            atan(y_coord(coord), x_coord(coord))
        )
    )
end

# Convenience functions for coordinate conversion
cartesian(coord::Coordinate{Frame, SphericalD{T, D}}) where {Frame, T, D} = spherical2cartesian(coord)
cartesian(coord::Coordinate{Frame, Spherical{T}}) where {Frame, T} = spherical2cartesian(coord)
cartesian(coord::Coordinate{Frame, CartesianD{T}}) where {Frame, T} = coord
cartesian(coord::Coordinate{Frame, Cartesian{T}}) where {Frame, T} = coord

spherical(coord::Coordinate{Frame, CartesianD{T}}) where {Frame, T} = cartesian2spherical(coord)
spherical(coord::Coordinate{Frame, Cartesian{T}}) where {Frame, T} = cartesian2spherical(coord)
spherical(coord::Coordinate{Frame, SphericalD{T, D}}) where {Frame, T, D} = coord
spherical(coord::Coordinate{Frame, Spherical{T}}) where {Frame, T} = coord

# Add distance to unit coordinates
function add_distance(d::D, coord::Coordinate{<:AbstractFrame, Spherical{T}}) where {D, T}
    return Coordinate(coord.frame, 
        SphericalD(
            coord.representation.latitude,
            coord.representation.longitude,
            d
        )
    )
end

# Distance already exists - replace it
function add_distance(d::D, coord::Coordinate{<:AbstractFrame, SphericalD{T, X}}) where {D, T, X}
    return Coordinate(coord.frame, 
        SphericalD(
            coord.representation.latitude,
            coord.representation.longitude,
            d
        )
    )
end

# Add distance to unit Cartesian coordinates
function add_distance(d::D, coord::Coordinate{<:AbstractFrame, Cartesian{T}}) where {D, T}
    return Coordinate(coord.frame, 
        CartesianD(
            SVector(
                x_coord(coord) * d,
                y_coord(coord) * d,
                z_coord(coord) * d
            )
        )
    )
end

# Replace distance in CartesianD coordinates
function add_distance(d::D, coord::Coordinate{<:AbstractFrame, CartesianD{T}}) where {D, T}
    previous_distance = dist(coord)
    return Coordinate(coord.frame, 
        CartesianD(
            SVector(
                x_coord(coord) * d / previous_distance,
                y_coord(coord) * d / previous_distance,
                z_coord(coord) * d / previous_distance
            )
        )
    )
end

# Remove distance from CartesianD coordinates (convert to unit vectors)
function remove_distance(coord::Coordinate{<:AbstractFrame, CartesianD{T}}) where {T}
    return Coordinate(coord.frame, 
        Cartesian(coord.representation.v)  # Cartesian constructor normalizes automatically
    )
end

# Remove distance from SphericalD coordinates
function remove_distance(coord::Coordinate{<:AbstractFrame, SphericalD{T, D}}) where {T, D}
    return Coordinate(coord.frame, 
        Spherical(
            lat(coord),
            lon(coord)
        )
    )
end

# No-op for coordinates that already don't have distance
remove_distance(coord::Coordinate{<:AbstractFrame, Cartesian{T}}) where {T} = coord
remove_distance(coord::Coordinate{<:AbstractFrame, Spherical{T}}) where {T} = coord

# Convenience function to create distance adder
function distance(d::D) where {D}
    coord -> return add_distance(d, coord)
end

function distance(d::Nothing)
    coord -> return remove_distance(coord)
end

