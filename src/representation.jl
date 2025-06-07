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
- `latitude::T`: Latitude angle
- `longitude::T`: Longitude angle

# Examples
```julia
s = Spherical(π/4, π/6)  # 45° latitude, 30° longitude
```
"""
struct Spherical{T} <: AbstractSphericalRepresentation
    latitude::T
    longitude::T
    function Spherical{T}(latitude, longitude) where {T}
        new{T}(latitude, longitude)
    end
end
Spherical(lat::T, lon::T) where {T} = Spherical{T}(lat, lon)
Spherical(lat, lon) = Spherical(promote(lat, lon)...)
Spherical{F}(c::T) where {F,T<:AbstractRepresentation} = convert(Spherical{F}, c)
"""
    dist(representation)

Get the distance component of a representation. Returns `1.0` for unit sphere representations.
"""
dist(representation::Spherical) = 1.0

"""
    SphericalD{T,D} <: AbstractSphericalRepresentation

Spherical coordinates with distance (latitude, longitude, distance).
Angles and distances can have different types and units.

# Fields
- `latitude::T`: Latitude angle
- `longitude::T`: Longitude angle  
- `distance::D`: Distance from origin

# Examples
```julia
s = SphericalD(π/4, π/6, 2.0)  # 45° latitude, 30° longitude, distance 2
# Could also be: SphericalD(12.3u"degree", 44u"degree", 23u"AU")
```
"""
struct SphericalD{T,D} <: AbstractSphericalRepresentation
    latitude::T
    longitude::T
    distance::D
    function SphericalD{T,D}(latitude, longitude, distance) where {T,D}
        new{T,D}(latitude, longitude, distance)
    end
end
SphericalD(lat::T, lon::T, dist::D) where {T, D} = SphericalD{T,D}(lat, lon, dist)
SphericalD(lat, lon, dist) = SphericalD(promote(lat, lon)..., dist)
SphericalD{T,D}(c::R, dist::Real) where {T,D,R<:AbstractRepresentation} = convert(SphericalD{T,D}, SphericalD(c, dist))
"""
    dist(representation)

Get the distance component of a representation. Returns `1.0` for unit sphere representations.
"""
dist(representation::SphericalD) = representation.distance
dist(representation::Spherical) = one(eltype(typeof(representation)))


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
dist(representation::Cartesian) = one(eltype(typeof(representation))) # Unit sphere distance


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
# Distance methods for cartesian representations
dist(representation::CartesianD) = sqrt(representation.x^2 + representation.y^2 + representation.z^2)



############################## COORDINATE CONVERSION HELPERS ##############################

# Helper functions to avoid code duplication
_cart_to_spherical_lat(x, y, z) = atan(z, sqrt(x^2 + y^2))
_cart_to_spherical_lon(x, y, z) = atan(y, x)
_spherical_to_cart_x(lat, lon, r=1) = r * cos(lat) * cos(lon)
_spherical_to_cart_y(lat, lon, r=1) = r * cos(lat) * sin(lon)
_spherical_to_cart_z(lat, lon, r=1) = r * sin(lat)

############################## CONVERSIONS ##############################

# Spherical conversions
Spherical(c::AbstractCartesianRepresentation) = Spherical(_cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)), 
                                                          _cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)))
Spherical(s::AbstractSphericalRepresentation) = Spherical(lat(s), lon(s))
Spherical(s::Spherical) = s

# SphericalD conversions
SphericalD(s::AbstractSphericalRepresentation, d) = SphericalD(lat(s), lon(s), d)
SphericalD(s::SphericalD) = s
SphericalD(c::Cartesian, d) = SphericalD(_cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)),
                                         _cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)), d)
SphericalD(c::CartesianD) = SphericalD(_cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c)),
                                       _cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c)),
                                       dist(c))

# Cartesian conversions  
Cartesian(s::AbstractSphericalRepresentation) = Cartesian(_spherical_to_cart_x(lat(s), lon(s)),
                                                          _spherical_to_cart_y(lat(s), lon(s)),
                                                          _spherical_to_cart_z(lat(s), lon(s)))
Cartesian(c::AbstractCartesianRepresentation) = Cartesian(x_coord(c), y_coord(c), z_coord(c))
Cartesian(c::Cartesian) = c

# CartesianD conversions
CartesianD(s::Spherical, d) = CartesianD(_spherical_to_cart_x(lat(s), lon(s), d),
                                         _spherical_to_cart_y(lat(s), lon(s), d),
                                         _spherical_to_cart_z(lat(s), lon(s), d))
CartesianD(s::SphericalD) = CartesianD(_spherical_to_cart_x(lat(s), lon(s), dist(s)),
                                       _spherical_to_cart_y(lat(s), lon(s), dist(s)),
                                       _spherical_to_cart_z(lat(s), lon(s), dist(s)))
CartesianD(c::Cartesian, d) = CartesianD(d * x_coord(c), d * y_coord(c), d * z_coord(c))
CartesianD(c::CartesianD) = c

############################## BASE.CONVERT METHODS ##############################

# Identity conversions
Base.convert(::Type{T}, c::T) where {T<:AbstractRepresentation} = c

# Spherical conversions
Base.convert(::Type{Spherical{T}}, c::AbstractCartesianRepresentation) where {T} = 
    Spherical{T}(T(_cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c))), 
                 T(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c))))

Base.convert(::Type{Spherical{T}}, s::AbstractSphericalRepresentation) where {T} = 
    Spherical{T}(T(lat(s)), T(lon(s)))

# SphericalD conversions  
Base.convert(::Type{SphericalD{T,D}}, s::AbstractSphericalRepresentation) where {T,D} =
    SphericalD{T,D}(T(lat(s)), T(lon(s)), D(hasfield(typeof(s), :distance) ? dist(s) : 1))

Base.convert(::Type{SphericalD{T,D}}, c::AbstractCartesianRepresentation) where {T,D} =
    SphericalD{T,D}(T(_cart_to_spherical_lat(x_coord(c), y_coord(c), z_coord(c))),
                    T(_cart_to_spherical_lon(x_coord(c), y_coord(c), z_coord(c))),
                    D(isa(c, CartesianD) ? dist(c) : 1))

# Cartesian conversions
Base.convert(::Type{Cartesian{T}}, s::AbstractSphericalRepresentation) where {T} =
    Cartesian{T}(T(_spherical_to_cart_x(lat(s), lon(s))),
                 T(_spherical_to_cart_y(lat(s), lon(s))),
                 T(_spherical_to_cart_z(lat(s), lon(s))))

Base.convert(::Type{Cartesian{T}}, c::AbstractCartesianRepresentation) where {T} =
    Cartesian{T}(T(x_coord(c)), T(y_coord(c)), T(z_coord(c)))

# CartesianD conversions
Base.convert(::Type{CartesianD{T}}, s::SphericalD) where {T} = 
    CartesianD{T}(T(_spherical_to_cart_x(lat(s), lon(s), dist(s))),
                  T(_spherical_to_cart_y(lat(s), lon(s), dist(s))),
                  T(_spherical_to_cart_z(lat(s), lon(s), dist(s))))

Base.convert(::Type{CartesianD{T}}, s::Spherical) where {T} =
    CartesianD{T}(T(_spherical_to_cart_x(lat(s), lon(s))),
                  T(_spherical_to_cart_y(lat(s), lon(s))),
                  T(_spherical_to_cart_z(lat(s), lon(s))))

Base.convert(::Type{CartesianD{T}}, c::AbstractCartesianRepresentation) where {T} = 
    CartesianD{T}(T(x_coord(c)), T(y_coord(c)), T(z_coord(c)))

############################## EQUALITY AND HASHING ##############################

"""
    ==(a::AbstractRepresentation, b::AbstractRepresentation)

Compare two coordinate representations for equality.
Representations are equal if they represent the same point in space,
regardless of coordinate system.

# Design Notes
- Same-type comparisons use field equality for efficiency
- Cross-type comparisons convert to canonical Cartesian form  
- Representations with/without distance are equal only if distance ≈ 1
- Hash functions ensure `isequal(x,y) ⟹ hash(x) == hash(y)` as required for Sets/Dicts

# Examples
```julia
s = Spherical(π/4, π/6)
c = Cartesian(s)
@assert s == c  # Different types, same point

sd = SphericalD(π/4, π/6, 1.0)  
@assert s == sd  # Unit distance matches unit sphere
```
"""
function Base.:(==)(a::T, b::T) where {T<:AbstractRepresentation}
    # Same type, compare fields directly
    return _fields_equal(a, b)
end

function Base.:(==)(a::AbstractRepresentation, b::AbstractRepresentation)
    # Different types - need to be careful about distance vs unit sphere
    _has_distance(::Union{SphericalD, CartesianD}) = true
    _has_distance(::AbstractRepresentation) = false
    
    # Both have distance or both are on unit sphere
    if _has_distance(a) == _has_distance(b)
        if _has_distance(a)
            # Both have distance, convert to CartesianD and compare
            ca, cb = CartesianD(a), CartesianD(b)
            return x_coord(ca) == x_coord(cb) && y_coord(ca) == y_coord(cb) && z_coord(ca) == z_coord(cb)
        else
            # Both on unit sphere, convert to Cartesian and compare
            ca, cb = Cartesian(a), Cartesian(b)
            return x_coord(ca) == x_coord(cb) && y_coord(ca) == y_coord(cb) && z_coord(ca) == z_coord(cb)
        end
    else
        # One has distance, one doesn't - they can only be equal if distance is 1
        if _has_distance(a)
            return dist(a) == 1 && Cartesian(a) == Cartesian(b)
        else
            return dist(b) == 1 && Cartesian(a) == Cartesian(b)
        end
    end
end

# Helper function to compare fields of same type
_fields_equal(a::Spherical, b::Spherical) = (a.latitude == b.latitude) && (a.longitude == b.longitude)
_fields_equal(a::SphericalD, b::SphericalD) = (a.latitude == b.latitude) && (a.longitude == b.longitude) && (a.distance == b.distance)
_fields_equal(a::Cartesian, b::Cartesian) = (a.x == b.x) && (a.y == b.y) && (a.z == b.z)
_fields_equal(a::CartesianD, b::CartesianD) = (a.x == b.x) && (a.y == b.y) && (a.z == b.z)

"""
    isequal(a::AbstractRepresentation, b::AbstractRepresentation)

Test whether two coordinate representations are equal, handling special cases like NaN.
"""
Base.isequal(a::T, b::T) where {T<:AbstractRepresentation} = _fields_isequal(a, b)

function Base.isequal(a::AbstractRepresentation, b::AbstractRepresentation)
    _has_distance(::Union{SphericalD, CartesianD}) = true
    _has_distance(::AbstractRepresentation) = false
    
    # Both have distance or both are on unit sphere
    if _has_distance(a) == _has_distance(b)
        if _has_distance(a)
            # Both have distance, convert to CartesianD and compare
            ca, cb = CartesianD(a), CartesianD(b)
            return isequal(x_coord(ca), x_coord(cb)) && isequal(y_coord(ca), y_coord(cb)) && isequal(z_coord(ca), z_coord(cb))
        else
            # Both on unit sphere, convert to Cartesian and compare
            ca, cb = Cartesian(a), Cartesian(b)
            return isequal(x_coord(ca), x_coord(cb)) && isequal(y_coord(ca), y_coord(cb)) && isequal(z_coord(ca), z_coord(cb))
        end
    else
        # One has distance, one doesn't - they can only be equal if distance is 1
        if _has_distance(a)
            return isequal(dist(a), 1) && isequal(Cartesian(a), Cartesian(b))
        else
            return isequal(dist(b), 1) && isequal(Cartesian(a), Cartesian(b))
        end
    end
end

# Helper function for isequal on same types
_fields_isequal(a::Spherical, b::Spherical) = isequal(a.latitude, b.latitude) && isequal(a.longitude, b.longitude)
_fields_isequal(a::SphericalD, b::SphericalD) = isequal(a.latitude, b.latitude) && isequal(a.longitude, b.longitude) && isequal(a.distance, b.distance)
_fields_isequal(a::Cartesian, b::Cartesian) = isequal(a.x, b.x) && isequal(a.y, b.y) && isequal(a.z, b.z)
_fields_isequal(a::CartesianD, b::CartesianD) = isequal(a.x, b.x) && isequal(a.y, b.y) && isequal(a.z, b.z)

"""
    hash(representation::AbstractRepresentation, h::UInt)

Compute hash code for coordinate representations.
Ensures that representations comparing equal have the same hash.
Uses Cartesian coordinates as canonical form for cross-type consistency.
"""
function Base.hash(rep::AbstractRepresentation, h::UInt)
    # Convert to Cartesian for canonical representation
    # This ensures that equal representations have equal hashes
    c = Cartesian(rep)
    h = hash(x_coord(c), h)
    h = hash(y_coord(c), h)
    h = hash(z_coord(c), h)
    return hash(:AbstractRepresentation, h)
end

# Specialized hash for types with distance to distinguish from unit sphere
function Base.hash(rep::SphericalD, h::UInt)
    c = CartesianD(rep)  # Preserve distance information
    h = hash(x_coord(c), h)
    h = hash(y_coord(c), h)
    h = hash(z_coord(c), h)
    return hash(:RepresentationWithDistance, h)
end

function Base.hash(rep::CartesianD, h::UInt)
    h = hash(rep.x, h)
    h = hash(rep.y, h)
    h = hash(rep.z, h)
    return hash(:RepresentationWithDistance, h)
end

############################## APPROXIMATE EQUALITY ##############################

"""
    isapprox(a::AbstractRepresentation, b::AbstractRepresentation; kwargs...)

Test whether two coordinate representations are approximately equal.
This is crucial for floating-point coordinate comparisons.

# Arguments
- `rtol::Real=Base.rtoldefault(...)`: relative tolerance
- `atol::Real=0`: absolute tolerance  
- `nans::Bool=false`: whether NaN values are considered equal

# Examples
```julia
s1 = Spherical(π/4, π/6)
s2 = Spherical(π/4 + 1e-15, π/6)
@assert isapprox(s1, s2)
```
"""
function Base.isapprox(a::T, b::T; kwargs...) where {T<:AbstractRepresentation}
    return _fields_isapprox(a, b; kwargs...)
end

function Base.isapprox(a::AbstractRepresentation, b::AbstractRepresentation; kwargs...)
    _has_distance(::Union{SphericalD, CartesianD}) = true
    _has_distance(::AbstractRepresentation) = false
    
    # Both have distance or both are on unit sphere
    if _has_distance(a) == _has_distance(b)
        if _has_distance(a)
            # Both have distance, convert to CartesianD and compare
            ca, cb = CartesianD(a), CartesianD(b)
            return isapprox(x_coord(ca), x_coord(cb); kwargs...) && 
                   isapprox(y_coord(ca), y_coord(cb); kwargs...) && 
                   isapprox(z_coord(ca), z_coord(cb); kwargs...)
        else
            # Both on unit sphere, convert to Cartesian and compare
            ca, cb = Cartesian(a), Cartesian(b)
            return isapprox(x_coord(ca), x_coord(cb); kwargs...) && 
                   isapprox(y_coord(ca), y_coord(cb); kwargs...) && 
                   isapprox(z_coord(ca), z_coord(cb); kwargs...)
        end
    else
        # One has distance, one doesn't - they can only be equal if distance ≈ 1
        if _has_distance(a)
            return isapprox(dist(a), 1; kwargs...) && isapprox(Cartesian(a), Cartesian(b); kwargs...)
        else
            return isapprox(dist(b), 1; kwargs...) && isapprox(Cartesian(a), Cartesian(b); kwargs...)
        end
    end
end

# Helper functions for isapprox on same types

# Generic angular distance that works with any number-like type
@inline function _angular_distance(a, b)
    T = promote_type(typeof(a), typeof(b))
    diff = abs(a - b)
    
    # Get 2π in the correct type/units
    # This works for Float64, Dual, Complex, Unitful, etc.
    period = 2 * convert(T, π)
    half_period = period / 2
    
    # Use the shorter arc distance
    return diff > half_period ? period - diff : diff
end

# Generic pole detection that works with any number type
@inline function _near_pole(lat, tolerance=nothing)
    T = typeof(lat)
    
    # Default tolerance based on type
    if tolerance === nothing
        tolerance = if T <: AbstractFloat
            sqrt(eps(T))  # Type-appropriate epsilon
        else
            sqrt(eps(float(real(T))))  # Handle Complex, Dual, etc.
        end
    end
    
    # Get π/2 in the correct type
    half_pi = convert(T, π) / 2
    
    return abs(abs(lat) - half_pi) < tolerance
end

function _fields_isapprox(a::Spherical, b::Spherical; kwargs...)
    # Fast latitude check first (most likely to fail)
    isapprox(a.latitude, b.latitude; kwargs...) || return false
    
    # Special case: at poles, longitude is undefined - always matches
    (_near_pole(a.latitude) || _near_pole(b.latitude)) && return true
    
    # Generic longitude wraparound comparison
    angular_dist = _angular_distance(a.longitude, b.longitude)
    
    # Get zero in the correct type for comparison
    zero_val = zero(typeof(angular_dist))
    return isapprox(angular_dist, zero_val; kwargs...)
end

function _fields_isapprox(a::SphericalD, b::SphericalD; kwargs...)
    # Fast checks with early return
    isapprox(a.latitude, b.latitude; kwargs...) || return false
    isapprox(a.distance, b.distance; kwargs...) || return false
    
    # Pole handling
    (_near_pole(a.latitude) || _near_pole(b.latitude)) && return true
    
    # Generic longitude comparison
    angular_dist = _angular_distance(a.longitude, b.longitude)
    zero_val = zero(typeof(angular_dist))
    return isapprox(angular_dist, zero_val; kwargs...)
end

function _fields_isapprox(a::Cartesian, b::Cartesian; kwargs...)
    return isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)
end

function _fields_isapprox(a::CartesianD, b::CartesianD; kwargs...)
    return isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)
end

############################## ADDITIONAL UTILITY FUNCTIONS ##############################

"""
    distance(from::AbstractRepresentation, to::AbstractRepresentation)

Calculate the angular distance between two coordinate representations on the unit sphere.
"""
function distance(from::AbstractRepresentation, to::AbstractRepresentation)
    c1 = Cartesian(from)
    c2 = Cartesian(to)
    dot_product = x_coord(c1) * x_coord(c2) + y_coord(c1) * y_coord(c2) + z_coord(c1) * z_coord(c2)
    return acos(clamp(dot_product, -1, 1))
end

"""
    norm(representation::AbstractCartesianRepresentation)

Calculate the Euclidean norm of a Cartesian representation.
"""
norm(representation::AbstractCartesianRepresentation) = sqrt(x_coord(representation)^2 + y_coord(representation)^2 + z_coord(representation)^2)

# Pretty printing
Base.show(io::IO, s::Spherical) = print(io, "Spherical(lat=$(s.latitude), lon=$(s.longitude))")
Base.show(io::IO, s::SphericalD) = print(io, "SphericalD{$(typeof(s.latitude)),$(typeof(s.distance))}(lat=$(s.latitude), lon=$(s.longitude), dist=$(s.distance))")
Base.show(io::IO, c::Cartesian) = print(io, "Cartesian(x=$(c.x), y=$(c.y), z=$(c.z))")
Base.show(io::IO, c::CartesianD) = print(io, "CartesianD(x=$(c.x), y=$(c.y), z=$(c.z))")










