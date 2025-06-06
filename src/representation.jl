abstract type AbstractRepresentation end 
abstract type AbstractSphericalRepresentation <: AbstractRepresentation end 
abstract type AbstractCartesianRepresentation <: AbstractRepresentation end  
Base.Broadcast.broadcastable(x::AbstractRepresentation) = Ref(x)



# Spherical coordinates without distance (unit sphere)
struct Spherical{T} <: AbstractSphericalRepresentation
    latitude::T
    longitude::T
    function Spherical(latitude::T, longitude::T) where {T}
        new{T}(latitude, longitude)
    end
end

# Spherical coordinates with distance
struct SphericalD{T, D} <: AbstractSphericalRepresentation
    latitude::T
    longitude::T
    distance::D
    function SphericalD(latitude::T, longitude::T, distance::D) where {T, D}
        new{T, D}(latitude, longitude, distance)
    end
end
SphericalD(s::AbstractSphericalRepresentation, d) = Spherical(s.latitude, s.longitude, d)

lon(representation::AbstractSphericalRepresentation) = representation.longitude
lat(representation::AbstractSphericalRepresentation) = representation.latitude
dist(representation::SphericalD{T}) where T = representation.distance

# Cartesian coordinates without explicit distance (unit vectors)
struct Cartesian{T} <: AbstractCartesianRepresentation
    x::T
    y::T
    z::T
    function Cartesian(x::T, y::T, z::T) where {T}
        n = sqrt(x^2 + y^2 + z^2)
        new{T}(x / n, y / n, z/  n)
    end
end

# Cartesian coordinates with distance
struct CartesianD{T} <: AbstractCartesianRepresentation
    x::T
    y::T
    z::T
    function CartesianD(x::T, y::T, z::T) where {T}
        new{T}(x,y,z)
    end
end

x_coord(representation::AbstractCartesianRepresentation) = representation.x
y_coord(representation::AbstractCartesianRepresentation) = representation.y
z_coord(representation::AbstractCartesianRepresentation) = representation.z
dist(representation::CartesianD{T}) where T = sqrt(representation.x^2 + representation.y^2 + representation.z^2)

function Cartesian(s::AbstractSphericalRepresentation)
    return Cartesian(
        cos(lat(coord)) * cos(lon(coord)),
        cos(lat(coord)) * sin(lon(coord)),
        sin(lat(coord))
    )
end

function CartesianD(s::SphericalD)
    return Cartesian(
        dist(s) * cos(lat(coord)) * cos(lon(coord)),
        dist(s) * cos(lat(coord)) * sin(lon(coord)),
        dist(s) * sin(lat(coord))
    )
end

function SphericalD(c::CartesianD)
    return SphericalD(
        atan(z_coord(coord), sqrt(x_coord(coord)^2 + y_coord(coord)^2)),
        atan(y_coord(coord), x_coord(coord)),
        sqrt(x_coord(coord)^2 + y_coord(coord)^2 + z_coord(coord)^2)
    )
end

function Spherical(c :: AbstractCartesianRepresentation)
    return SphericalD(
        atan(z_coord(coord), sqrt(x_coord(coord)^2 + y_coord(coord)^2)),
        atan(y_coord(coord), x_coord(coord))
    )
end
