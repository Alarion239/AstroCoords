"""
    AbstractRepresentation

Abstract base type for all coordinate representations in astronomical coordinate systems.
All concrete representations should inherit from this type either directly or through
intermediate abstract types.
"""
abstract type AbstractRepresentation end 
Base.Broadcast.broadcastable(x::AbstractRepresentation) = Ref(x)

############################## SPHERICAL REPRESENTATIONS ##############################


"""
    Spherical{T, D} <: AbstractRepresentation

Spherical coordinates (latitude, longitude, distance).

# Fields
- `longitude::T`: Longitude angle
- `latitude::T`: Latitude angle
- `distance::D`: Distance from origin
# Examples
```julia
s = Spherical(π/6, π/4)  # 30° longitude, 45° latitude
s = Spherical(π/6, π/4, 1.0)  # 30° longitude, 45° latitude, distance 1.0
```
"""
struct Spherical{T, D} <: AbstractRepresentation
    longitude::T
    latitude::T
    distance::D
    function Spherical{T,D}(longitude, latitude, distance) where {T,D}
        lat_norm, lon_norm = _normalize_lat_lon(latitude, longitude)
        new{T,D}(lon_norm, lat_norm, distance)
    end
end
Spherical(lon::T, lat::T) where {T} = Spherical{T, Float64}(lon, lat, 1.0)
Spherical(lon, lat) = Spherical(promote(lon, lat)..., 1.0)
Spherical{F}(c::AbstractRepresentation) where {F} = convert(Spherical{F, Float64}, c)
Spherical(lon::T, lat::T, dist::D) where {T, D} = Spherical{T,D}(lon, lat, dist)
Spherical(lon, lat, dist) = Spherical(promote(lon, lat)..., dist)
Spherical{T,D}(c::R, dist::Real) where {T,D,R<:AbstractRepresentation} = convert(Spherical{T,D}, Spherical(c, dist))

# Optimized helper functions for angle normalization
function _normalize_lat_lon(lat::T, lon::T) where T
    # Use rem2pi for better precision than mod(x, 2π)
    lat_normalized = rem2pi(lat, RoundDown)  # Normalize to [0, 2π] with higher precision
    lon_normalized = rem2pi(lon, RoundNearest)  # Normalize to [-π, π] directly
    
    # Early return if latitude is already in valid range
    if lat_normalized <= π/2
        return lat_normalized, lon_normalized
    elseif lat_normalized >= 3π/2
        return lat_normalized - 2π, lon_normalized
    else
        # Handle pole crossing (π/2 < lat < 3π/2)
        lat_norm = π - lat_normalized
        # Use muladd for better precision: lon + π
        lon_norm = rem2pi(muladd(lon_normalized, 1, π), RoundNearest)
        return lat_norm, lon_norm
    end
end

"""
    lon(representation::Spherical)

Get the longitude component of a spherical representation.
"""
lon(representation::Spherical) = representation.longitude

"""
    lat(representation::Spherical)

Get the latitude component of a spherical representation.
"""
lat(representation::Spherical) = representation.latitude





############################## CARTESIAN REPRESENTATIONS ##############################


"""
    Cartesian{T} <: AbstractRepresentation

Normalized Cartesian coordinates on a unit sphere (x, y, z with ||r|| = 1).

# Fields
- `x::T`: X-coordinate (normalized)
- `y::T`: Y-coordinate (normalized)
- `z::T`: Z-coordinate (normalized)

# Examples
```julia
c = Cartesian(1.0, 1.0, 1.0)  # Automatically normalized to unit length
```
"""
struct Cartesian{T} <: AbstractRepresentation
    x::T
    y::T
    z::T
    distance::T
    function Cartesian{T}(x, y, z, distance) where {T}
        new{T}(x, y, z, distance)
    end
end
Cartesian(x::D, y::D, z::D, distance::D) where {D} = Cartesian{D}(x, y, z, distance)  
Cartesian(x::D, y::D, z::D, distance::D) where {D <: Real} = Cartesian{float(D)}(x, y, z, distance)  
Cartesian(x, y, z, distance) = Cartesian(promote(x, y, z, distance)...)
Cartesian{F}(c::T) where {F,T<:AbstractRepresentation} = convert(Cartesian{F, typeof(dist(c))}, c)
Cartesian{F}(c::Cartesian{T}) where {F,T} = Cartesian{F}(c.x, c.y, c.z, c.distance)
Cartesian{F}(c::Cartesian{T}, distance::D) where {F,T,D} = Cartesian{F}(c.x, c.y, c.z, distance)

function Cartesian(x, y, z) 
    distance = sqrt(x^2 + y^2 + z^2)
    if iszero(distance)
        throw(ArgumentError("Cannot normalize zero vector"))
    end
    Cartesian(x / distance, y / distance, z / distance, distance)
end

"""
    x_coord(representation::Cartesian)

Get the x-coordinate component of a cartesian representation.
"""
x_coord(representation::Cartesian) = representation.x

"""
    y_coord(representation::Cartesian)

Get the y-coordinate component of a cartesian representation.
"""
y_coord(representation::Cartesian) = representation.y

"""
    z_coord(representation::Cartesian)

Get the z-coordinate component of a cartesian representation.
"""
z_coord(representation::Cartesian) = representation.z

"""
    dist(representation::AbstractRepresentation)

Get the distance component of a coordinate representation.

# Returns
- For `Spherical{T,D}`: Returns the stored distance value of type `D`
- For `Cartesian{D}`: Returns the stored distance value of type `D`

# Examples
```julia
julia> dist(Spherical(π/4, π/6))
1.0

julia> dist(Cartesian(3.0, 4.0, 0.0))
5.0
```
"""

dist(representation::Spherical) = representation.distance
dist(representation::Cartesian) = representation.distance



############################## COORDINATE CONVERSION HELPERS ##############################

# Helper functions to avoid code duplication
_cart_to_spherical_lon(x, y, z) = atan(y, x)
_cart_to_spherical_lat(x, y, z) = atan(z, sqrt(x^2 + y^2))
_spherical_to_cart_x(lon, lat, r) = r * cos(lat) * cos(lon)
_spherical_to_cart_y(lon, lat, r) = r * cos(lat) * sin(lon)
_spherical_to_cart_z(lon, lat, r) = r * sin(lat)

############################## CONVERSIONS ##############################

# Spherical conversions
Spherical(c::Cartesian) = Spherical(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)), 
                                    _cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)),
                                    dist(c))
Spherical(s::Spherical) = s

# Cartesian conversions  
Cartesian(s::Spherical) = Cartesian(_spherical_to_cart_x(lon(s), lat(s), dist(s)),
                                    _spherical_to_cart_y(lon(s), lat(s), dist(s)),
                                    _spherical_to_cart_z(lon(s), lat(s), dist(s)))
Cartesian(c::Cartesian) = c

############################## BASE.CONVERT METHODS ##############################

# Identity conversions
Base.convert(::Type{T}, c::T) where {T<:AbstractRepresentation} = c

# Spherical conversions
Base.convert(::Type{Spherical{T}}, c::Cartesian) where {T} = Spherical(c)
Base.convert(::Type{Spherical{T}}, s::Spherical) where {T} = Spherical(s)

# Cartesian conversions
Base.convert(::Type{Cartesian{T}}, s::Spherical) where {T} = Cartesian(s)
Base.convert(::Type{Cartesian{T}}, c::Cartesian) where {T} = Cartesian(c)

# Pretty printing
Base.show(io::IO, s::Spherical) = print(io, "Spherical{$(typeof(s.longitude))}(lon=$(s.longitude), lat=$(s.latitude))")
Base.show(io::IO, c::Cartesian) = print(io, "Cartesian{$(typeof(c.x))}(x=$(c.x), y=$(c.y), z=$(c.z))")











