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
    AbstractSphericalRepresentation <: AbstractRepresentation

Abstract type for spherical coordinate representations (latitude, longitude, [distance]).
"""
abstract type AbstractSphericalRepresentation <: AbstractRepresentation end 

"""
    lon(representation::AbstractSphericalRepresentation)

Get the longitude component of a spherical representation.
"""
lon(representation::AbstractSphericalRepresentation) = representation.longitude

"""
    lat(representation::AbstractSphericalRepresentation)

Get the latitude component of a spherical representation.
"""
lat(representation::AbstractSphericalRepresentation) = representation.latitude



"""
    Spherical{T} <: AbstractSphericalRepresentation

Spherical coordinates on a unit sphere (latitude, longitude).

# Fields
- `longitude::T`: Longitude angle
- `latitude::T`: Latitude angle

# Examples
```julia
s = Spherical(π/6, π/4)  # 30° longitude, 45° latitude
```
"""
struct Spherical{T} <: AbstractSphericalRepresentation
    longitude::T
    latitude::T
    function Spherical{T}(longitude, latitude) where {T}
        new{T}(longitude, latitude)
    end
end
Spherical(lon::T, lat::T) where {T} = Spherical{T}(lon, lat)
Spherical(lon, lat) = Spherical(promote(lon, lat)...)
Spherical{F}(c::AbstractRepresentation) where {F} = convert(Spherical{F}, c)



"""
    SphericalD{T,D} <: AbstractSphericalRepresentation

Spherical coordinates with distance (latitude, longitude, distance).
Angles and distances can have different types and units.

# Fields
- `longitude::T`: Longitude angle
- `latitude::T`: Latitude angle
- `distance::D`: Distance from origin

# Examples
```julia
s = SphericalD(π/6, π/4, 2.0)  # 30° longitude, 45° latitude, distance 2
# Could also be: SphericalD(12.3u"degree", 44u"degree", 23u"AU")
```
"""
struct SphericalD{T,D} <: AbstractSphericalRepresentation
    longitude::T
    latitude::T
    distance::D
    function SphericalD{T,D}(longitude, latitude, distance) where {T,D}
        new{T,D}(longitude, latitude, distance)
    end
end
SphericalD(lon::T, lat::T, dist::D) where {T, D} = SphericalD{T,D}(lon, lat, dist)
SphericalD(lon, lat, dist) = SphericalD(promote(lon, lat)..., dist)
SphericalD{T,D}(c::R, dist::Real) where {T,D,R<:AbstractRepresentation} = convert(SphericalD{T,D}, SphericalD(c, dist))



############################## CARTESIAN REPRESENTATIONS ##############################

"""
    AbstractCartesianRepresentation <: AbstractRepresentation

Abstract type for Cartesian coordinate representations (x, y, z, [with normalization]).
"""
abstract type AbstractCartesianRepresentation <: AbstractRepresentation end  

"""
    x_coord(representation::AbstractCartesianRepresentation)

Get the x-coordinate component of a cartesian representation.
"""
x_coord(representation::AbstractCartesianRepresentation) = representation.x

"""
    y_coord(representation::AbstractCartesianRepresentation)

Get the y-coordinate component of a cartesian representation.
"""
y_coord(representation::AbstractCartesianRepresentation) = representation.y

"""
    z_coord(representation::AbstractCartesianRepresentation)

Get the z-coordinate component of a cartesian representation.
"""
z_coord(representation::AbstractCartesianRepresentation) = representation.z

"""
    Cartesian{T} <: AbstractCartesianRepresentation

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
struct Cartesian{T} <: AbstractCartesianRepresentation
    x::T
    y::T
    z::T
    function Cartesian{T}(x, y, z) where {T}
        n = sqrt(x^2 + y^2 + z^2)
        if iszero(n)
            throw(ArgumentError("Cannot normalize zero vector"))
        end
        new{T}(x / n, y / n, z / n)
    end
end
Cartesian(x::T, y::T, z::T) where {T} = Cartesian{T}(x, y, z)  
Cartesian(x::T, y::T, z::T) where {T <: Real} = Cartesian{float(T)}(x, y, z)  
Cartesian(x, y, z) = Cartesian(promote(x, y, z)...)
Cartesian{F}(c::T) where {F,T<:AbstractRepresentation} = convert(Cartesian{F}, c)



"""
    CartesianD{T} <: AbstractCartesianRepresentation

Cartesian coordinates with explicit distance (x, y, z).

# Fields
- `x::T`: X-coordinate
- `y::T`: Y-coordinate
- `z::T`: Z-coordinate

# Examples
```julia
c = CartesianD(1.0, 2.0, 3.0)  # Cartesian coordinates with distance sqrt(14)
```
"""
struct CartesianD{T} <: AbstractCartesianRepresentation
    x::T
    y::T
    z::T
    function CartesianD{T}(x, y, z) where {T}
        new{T}(x,y,z)
    end
end
CartesianD(x::T, y::T, z::T) where {T} = CartesianD{T}(x, y, z)
CartesianD(x, y, z) = CartesianD(promote(x, y, z)...)
CartesianD{F}(c::T) where {F,T<:AbstractRepresentation} = convert(CartesianD{F}, c)

"""
    dist(representation::AbstractRepresentation)

Get the distance component of a coordinate representation.

# Returns
- For `Spherical{T}`: Returns `1.0` of type `Float64` (unit sphere)
- For `Cartesian{T}`: Returns unit length of the type `T` (unit sphere)
- For `SphericalD{T,D}`: Returns the stored distance value of type `D`
- For `CartesianD{T}`: Returns the computed Euclidean distance √(x² + y² + z²) of type `T`

# Examples
```julia
julia> dist(Spherical(π/4, π/6))
1.0

julia> dist(SphericalD(π/4, π/6, 5.0))
5.0

julia> dist(CartesianD(3.0, 4.0, 0.0))
5.0
```
"""
dist(representation::SphericalD) = representation.distance
dist(representation::CartesianD) = sqrt(representation.x^2 + representation.y^2 + representation.z^2)
dist(::Cartesian{T}) where {T} = one(T)
dist(representation::Spherical) = 1.0



############################## COORDINATE CONVERSION HELPERS ##############################

# Helper functions to avoid code duplication
_cart_to_spherical_lon(x, y, z) = atan(y, x)
_cart_to_spherical_lat(x, y, z) = atan(z, sqrt(x^2 + y^2))
_spherical_to_cart_x(lon, lat, r=1) = r * cos(lat) * cos(lon)
_spherical_to_cart_y(lon, lat, r=1) = r * cos(lat) * sin(lon)
_spherical_to_cart_z(lon, lat, r=1) = r * sin(lat)

############################## CONVERSIONS ##############################

# Spherical conversions
Spherical(c::AbstractCartesianRepresentation) = Spherical(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)), 
                                                          _cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)))
Spherical(s::AbstractSphericalRepresentation) = Spherical(lon(s), lat(s))
Spherical(s::Spherical) = s

# SphericalD conversions
SphericalD(s::AbstractSphericalRepresentation; distance = dist(s)) = SphericalD(lon(s), lat(s), distance)
SphericalD(s::SphericalD; distance = dist(s)) = SphericalD(lon(s), lat(s), distance)
SphericalD(c::Cartesian{T}; distance = one(T)) where {T} = SphericalD(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)),
                                         _cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)), distance)
SphericalD(c::CartesianD; distance = dist(c)) = SphericalD(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)),
                                       _cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)),
                                       distance)
SphericalD(d) = (x::AbstractRepresentation) -> SphericalD(x, distance = d)

# Cartesian conversions  
Cartesian(s::AbstractSphericalRepresentation) = Cartesian(_spherical_to_cart_x(lon(s), lat(s)),
                                                          _spherical_to_cart_y(lon(s), lat(s)),
                                                          _spherical_to_cart_z(lon(s), lat(s)))
Cartesian(c::AbstractCartesianRepresentation) = Cartesian(x_coord(c), y_coord(c), z_coord(c))
Cartesian(c::Cartesian) = c

# CartesianD conversions
CartesianD(s::Spherical; distance = dist(s)) = CartesianD(_spherical_to_cart_x(lon(s), lat(s), distance),
                                         _spherical_to_cart_y(lon(s), lat(s), distance),
                                         _spherical_to_cart_z(lon(s), lat(s), distance))
CartesianD(s::SphericalD; distance = dist(s)) = CartesianD(_spherical_to_cart_x(lon(s), lat(s), distance),
                                       _spherical_to_cart_y(lon(s), lat(s), dist(s)),
                                       _spherical_to_cart_z(lon(s), lat(s), dist(s)))
CartesianD(c::Cartesian; distance = dist(c)) = CartesianD(distance * x_coord(c), distance * y_coord(c), distance * z_coord(c))
CartesianD(c::CartesianD; distance = dist(c)) = CartesianD(distance * x_coord(c), distance * y_coord(c), distance * z_coord(c))
CartesianD(d) = (x::AbstractRepresentation) -> CartesianD(x, distance = d)

############################## BASE.CONVERT METHODS ##############################

# Identity conversions
Base.convert(::Type{T}, c::T) where {T<:AbstractRepresentation} = c

# Spherical conversions
Base.convert(::Type{Spherical{T}}, c::AbstractCartesianRepresentation) where {T} = Spherical(c)
Base.convert(::Type{Spherical{T}}, s::AbstractSphericalRepresentation) where {T} = Spherical(s)
Base.convert(::Type{Spherical{T}}, s::Spherical) where {T} = Spherical(s)

# SphericalD conversions  
Base.convert(::Type{SphericalD{T,D}}, s::AbstractSphericalRepresentation) where {T,D} = SphericalD(s)
Base.convert(::Type{SphericalD{T,D}}, c::AbstractCartesianRepresentation) where {T,D} = SphericalD(c)

# Cartesian conversions
Base.convert(::Type{Cartesian{T}}, s::AbstractSphericalRepresentation) where {T} = Cartesian(s)
Base.convert(::Type{Cartesian{T}}, c::AbstractCartesianRepresentation) where {T} = Cartesian(c)

# CartesianD conversions
Base.convert(::Type{CartesianD{T}}, s::SphericalD) where {T} = CartesianD(s)
Base.convert(::Type{CartesianD{T}}, s::Spherical) where {T} = CartesianD(s)
Base.convert(::Type{CartesianD{T}}, c::AbstractCartesianRepresentation) where {T} = CartesianD(c)

# Pretty printing
Base.show(io::IO, s::Spherical) = print(io, "Spherical{$(typeof(s.longitude))}(lon=$(s.longitude), lat=$(s.latitude))")
Base.show(io::IO, s::SphericalD) = print(io, "SphericalD{$(typeof(s.longitude)),$(typeof(s.distance))}(lon=$(s.longitude), lat=$(s.latitude), dist=$(s.distance))")
Base.show(io::IO, c::Cartesian) = print(io, "Cartesian{$(typeof(c.x))}(x=$(c.x), y=$(c.y), z=$(c.z))")
Base.show(io::IO, c::CartesianD) = print(io, "CartesianD{$(typeof(c.x))}(x=$(c.x), y=$(c.y), z=$(c.z))")










