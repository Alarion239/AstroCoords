"""
    CoordinateAccessors

This module provides both generic and frame-specific coordinate component accessors for astronomical coordinate systems.
It follows the astropy convention but with Julia-optimized performance through direct representation access.

## Generic Coordinate Accessors
Generic accessors work with any coordinate regardless of frame type:
- `lon`, `lat`, `dist` - For spherical representations
- `x_coord`, `y_coord`, `z_coord`, `dist` - For Cartesian representations

## Frame-Specific Accessors
The accessors are specialized for specific frame types to provide natural naming conventions:
- ICRS: `ra`, `dec` (right ascension, declination)
- Galactic: `l`, `b` (galactic longitude, latitude)
- AltAz: `alt`, `az` (altitude, azimuth)
- HADec: `ha`, `dec` (hour angle, declination)
- FK4/FK5: `ra`, `dec` (right ascension, declination)
- Ecliptic: `λ`, `β` or `eclon`, `eclat` (ecliptic longitude, latitude)
- Supergalactic: `sgl`, `sgb` (supergalactic longitude, latitude)

All accessors are optimized for maximum performance using direct representation access and compile-time dispatch.
"""
#####
using AstroCoords: Coordinate, AbstractFrame, AbstractRepresentation, 
ICRS, Galactic, AltAz, FK4, FK5, FK4NoETerms, HADec, Ecliptic, Supergalactic, TEME, CIRS, ITRS, HCRS

# =============================================================================
# Generic Coordinate Component Accessors
# =============================================================================

"""
    lon(coordinate::Coordinate{AbstractFrame, Spherical})

Extract the longitude component from a coordinate with spherical representation.

# Arguments
- `coordinate`: Coordinate with spherical representation

# Returns
- Longitude component from the underlying representation
"""
lon(coordinate::Coordinate{AbstractFrame, Spherical}) = lon(coordinate.representation)

"""
    lat(coordinate::Coordinate{AbstractFrame, Spherical})

Extract the latitude component from a coordinate with spherical representation.

# Arguments
- `coordinate`: Coordinate with spherical representation

# Returns
- Latitude component from the underlying representation
"""
lat(coordinate::Coordinate{AbstractFrame, Spherical}) = lat(coordinate.representation)

"""
    dist(coordinate::Coordinate{AbstractFrame, AbstractRepresentation})

Extract the distance component from a coordinate with spherical representation that includes distance.

# Arguments
- `coordinate`: Coordinate with `AbstractRepresentation` representation

# Returns
- Distance component from the underlying representation
"""
dist(coordinate::Coordinate{AbstractFrame, AbstractRepresentation}) = dist(coordinate.representation)

"""
    x_coord(coordinate::Coordinate{AbstractFrame, Cartesian})

Extract the x-coordinate component from a coordinate with Cartesian representation.

# Arguments
- `coordinate`: Coordinate with Cartesian representation

# Returns
- X-component from the underlying representation
"""
x_coord(coordinate::Coordinate{AbstractFrame, Cartesian}) = x_coord(coordinate.representation)

"""
    y_coord(coordinate::Coordinate{AbstractFrame, Cartesian})

Extract the y-coordinate component from a coordinate with Cartesian representation.

# Arguments
- `coordinate`: Coordinate with Cartesian representation

# Returns
- Y-component from the underlying representation
"""
y_coord(coordinate::Coordinate{AbstractFrame, Cartesian}) = y_coord(coordinate.representation)

"""
    z_coord(coordinate::Coordinate{AbstractFrame, Cartesian})

Extract the z-coordinate component from a coordinate with Cartesian representation.

# Arguments
- `coordinate`: Coordinate with Cartesian representation

# Returns
- Z-component from the underlying representation
"""
z_coord(coordinate::Coordinate{AbstractFrame, Cartesian}) = z_coord(coordinate.representation)


# =============================================================================
# ICRS Frame Accessors (Right Ascension, Declination)
# =============================================================================

"""
    ra(c::Coordinate{F, R}) where {F::ICRS,R::Spherical}

Extract the right ascension from an ICRS coordinate.

# Arguments
- `c`: Coordinate in ICRS frame with spherical representation

# Returns
- Right ascension component (same type as coordinate representation)

# Examples
```julia-repl
julia> coord = Coordinate(ICRS(), SphericalD(deg2rad(30), deg2rad(45), 1.0))
julia> ra(coord)
0.5235987755982988
```
"""
ra(c::Coordinate{ICRS, Spherical}) = lon(c.representation)

"""
    dec(c::Coordinate{ICRS, Spherical})

Extract the declination from an ICRS coordinate.

# Arguments
- `c`: Coordinate in ICRS frame with spherical representation

# Returns
- Declination component (same type as coordinate representation)
"""
dec(c::Coordinate{ICRS, Spherical}) = lat(c.representation)

# =============================================================================
# Galactic Frame Accessors (Galactic Longitude, Galactic Latitude)
# =============================================================================

"""
    l(c::Coordinate{Galactic, Spherical})

Extract the galactic longitude from a Galactic coordinate.

# Arguments
- `c`: Coordinate in Galactic frame with spherical representation

# Returns
- Galactic longitude component (same type as coordinate representation)
"""
gal(c::Coordinate{Galactic, Spherical}) = lon(c.representation)

"""
    b(c::Coordinate{Galactic, Spherical})

Extract the galactic latitude from a Galactic coordinate.

# Arguments
- `c`: Coordinate in Galactic frame with spherical representation

# Returns
- Galactic latitude component (same type as coordinate representation)
"""
gab(c::Coordinate{Galactic, Spherical}) = lat(c.representation)

# =============================================================================
# AltAz Frame Accessors (Altitude, Azimuth)
# =============================================================================

"""
    alt(c::Coordinate{<:AltAz, Spherical})

Extract the altitude (elevation) from an AltAz coordinate.

# Arguments
- `c`: Coordinate in AltAz frame with spherical representation

# Returns
- Altitude component (same type as coordinate representation)
"""
alt(c::Coordinate{AltAz, Spherical}) = lat(c.representation)

"""
    az(c::Coordinate{<:AltAz, Spherical})

Extract the azimuth from an AltAz coordinate.

# Arguments
- `c`: Coordinate in AltAz frame with spherical representation

# Returns
- Azimuth component (same type as coordinate representation)
"""
az(c::Coordinate{AltAz, Spherical}) = lon(c.representation)

# =============================================================================
# HADec Frame Accessors (Hour Angle, Declination)
# =============================================================================

"""
    ha(c::Coordinate{<:HADec, Spherical})

Extract the hour angle from a HADec coordinate.

# Arguments
- `c`: Coordinate in HADec frame with spherical representation

# Returns
- Hour angle component (same type as coordinate representation)
"""
ha(c::Coordinate{HADec, Spherical}) = lon(c.representation)

"""
    dec(c::Coordinate{<:HADec, Spherical})

Extract the declination from a HADec coordinate.

# Arguments
- `c`: Coordinate in HADec frame with spherical representation

# Returns
- Declination component (same type as coordinate representation)
"""
dec(c::Coordinate{HADec, Spherical}) = lat(c.representation)

# =============================================================================
# FK4/FK5 Frame Accessors (Right Ascension, Declination)
# =============================================================================

"""
    ra(c::Coordinate{<:FK4, Spherical})

Extract the right ascension from an FK4 coordinate.

# Arguments
- `c`: Coordinate in FK4 frame with spherical representation

# Returns
- Right ascension component (same type as coordinate representation)
"""
ra(c::Coordinate{FK4, Spherical}) = lon(c.representation)

"""
    dec(c::Coordinate{<:FK4, Spherical})

Extract the declination from an FK4 coordinate.

# Arguments
- `c`: Coordinate in FK4 frame with spherical representation

# Returns
- Declination component (same type as coordinate representation)
"""
dec(c::Coordinate{FK4, Spherical}) = lat(c.representation)

"""
    ra(c::Coordinate{<:FK5, Spherical})

Extract the right ascension from an FK5 coordinate.

# Arguments
- `c`: Coordinate in FK5 frame with spherical representation

# Returns
- Right ascension component (same type as coordinate representation)
"""
ra(c::Coordinate{FK5, Spherical}) = lon(c.representation)

"""
    dec(c::Coordinate{<:FK5, Spherical})

Extract the declination from an FK5 coordinate.

# Arguments
- `c`: Coordinate in FK5 frame with spherical representation

# Returns
- Declination component (same type as coordinate representation)
"""
dec(c::Coordinate{FK5, Spherical}) = lat(c.representation)

"""
    ra(c::Coordinate{<:FK4NoETerms, Spherical})

Extract the right ascension from an FK4NoETerms coordinate.

# Arguments
- `c`: Coordinate in FK4NoETerms frame with spherical representation

# Returns
- Right ascension component (same type as coordinate representation)
"""
ra(c::Coordinate{FK4NoETerms, Spherical}) = lon(c.representation)

"""
    dec(c::Coordinate{<:FK4NoETerms, Spherical})

Extract the declination from an FK4NoETerms coordinate.

# Arguments
- `c`: Coordinate in FK4NoETerms frame with spherical representation

# Returns
- Declination component (same type as coordinate representation)
"""
dec(c::Coordinate{FK4NoETerms, Spherical}) = lat(c.representation)

# =============================================================================
# Ecliptic Frame Accessors (Ecliptic Longitude, Ecliptic Latitude)
# =============================================================================

"""
    λ(c::Coordinate{<:Ecliptic, Spherical})

Extract the ecliptic longitude (lambda) from an Ecliptic coordinate.

# Arguments
- `c`: Coordinate in Ecliptic frame with spherical representation

# Returns
- Ecliptic longitude component (same type as coordinate representation)
"""
λ(c::Coordinate{Ecliptic, Spherical}) = lon(c.representation)

"""
    β(c::Coordinate{<:Ecliptic, Spherical})

Extract the ecliptic latitude (beta) from an Ecliptic coordinate.

# Arguments
- `c`: Coordinate in Ecliptic frame with spherical representation

# Returns
- Ecliptic latitude component (same type as coordinate representation)
"""
β(c::Coordinate{Ecliptic, Spherical}) = lat(c.representation)

"""
    eclon(c::Coordinate{<:Ecliptic, Spherical})

Extract the ecliptic longitude from an Ecliptic coordinate (ASCII alternative to λ).

# Arguments
- `c`: Coordinate in Ecliptic frame with spherical representation

# Returns
- Ecliptic longitude component (same type as coordinate representation)
"""
eclon(c::Coordinate{Ecliptic, Spherical}) = lon(c.representation)

"""
    eclat(c::Coordinate{<:Ecliptic, Spherical})

Extract the ecliptic latitude from an Ecliptic coordinate (ASCII alternative to β).

# Arguments
- `c`: Coordinate in Ecliptic frame with spherical representation

# Returns
- Ecliptic latitude component (same type as coordinate representation)
"""
eclat(c::Coordinate{Ecliptic, Spherical}) = lat(c.representation)

# =============================================================================
# Supergalactic Frame Accessors
# =============================================================================

"""
    sgl(c::Coordinate{Supergalactic, Spherical})

Extract the supergalactic longitude from a Supergalactic coordinate.

# Arguments
- `c`: Coordinate in Supergalactic frame with spherical representation

# Returns
- Supergalactic longitude component (same type as coordinate representation)
"""
sgl(c::Coordinate{Supergalactic, Spherical}) = lon(c.representation)

"""
    sgb(c::Coordinate{Supergalactic, Spherical})

Extract the supergalactic latitude from a Supergalactic coordinate.

# Arguments
- `c`: Coordinate in Supergalactic frame with spherical representation

# Returns
- Supergalactic latitude component (same type as coordinate representation)
"""
sgb(c::Coordinate{Supergalactic, Spherical}) = lat(c.representation)

# =============================================================================
# High-Performance Property Access Interface
# =============================================================================

"""
    Base.getproperty(c::Coordinate, s::Symbol)

Provide property-style access to coordinate components and frame-specific accessors.

This method enables both generic coordinate access (`.lon`, `.lat`, `.dist`, `.x`, `.y`, `.z`)
and frame-specific access (`.ra`, `.dec`, `.l`, `.b`, etc.) through Julia's property syntax.

Performance is optimized using constant folding with simple if-else chains rather than
dictionaries or other dynamic dispatch mechanisms.

# Arguments
- `c`: Coordinate object
- `s`: Property symbol to access

# Returns
- Requested coordinate component

# Supported Properties

## Core Properties
- `:frame`: The coordinate frame
- `:representation`: The coordinate representation

## Generic Coordinate Properties
- `:lon`, `:lat`, `:dist`: Spherical coordinates
- `:x`, `:y`, `:z`: Cartesian coordinates (aliases for x_coord, y_coord, z_coord)
- `:x_coord`, `:y_coord`, `:z_coord`: Cartesian coordinates

## Frame-Specific Properties
- `:ra`, `:dec`: Right ascension, declination (ICRS, FK4, FK5, FK4NoETerms, HADec)
- `:l`, `:b`: Galactic longitude, latitude (Galactic)
- `:alt`, `:az`: Altitude, azimuth (AltAz)
- `:ha`: Hour angle (HADec)
- `:λ`, `:β`: Ecliptic longitude, latitude (Ecliptic, Unicode)
- `:eclon`, `:eclat`: Ecliptic longitude, latitude (Ecliptic, ASCII)
- `:sgl`, `:sgb`: Supergalactic longitude, latitude (Supergalactic)

# Examples
```julia-repl
julia> coord = Coordinate(ICRS(), SphericalD(deg2rad(30), deg2rad(45), 1.0))
julia> coord.ra     # Right ascension
julia> coord.dec    # Declination
julia> coord.lon    # Generic longitude
julia> coord.lat    # Generic latitude
```

# Throws
- `ArgumentError`: If the requested property is not available for this coordinate
"""
@inline function Base.getproperty(c::Coordinate, s::Symbol)
    # Core properties - direct field access
    if s === :frame
        return getfield(c, :frame)
    elseif s === :representation
        return getfield(c, :representation)
    
    # Generic coordinate properties (representation-dependent)
    elseif s === :lon
        return lon(c)
    elseif s === :lat
        return lat(c)
    elseif s === :dist
        return dist(c)  # Will throw for representations without distance
    elseif s === :x
        return x_coord(c)
    elseif s === :y
        return y_coord(c)
    elseif s === :z
        return z_coord(c)
    elseif s === :x_coord
        return x_coord(c)
    elseif s === :y_coord
        return y_coord(c)
    elseif s === :z_coord
        return z_coord(c)
    
    # Frame-specific properties (frame-dependent)
    elseif s === :ra
        return ra(c)
    elseif s === :dec
        return dec(c)
    elseif s === :l
        return l(c)
    elseif s === :b
        return b(c)
    elseif s === :alt
        return alt(c)
    elseif s === :az
        return az(c)
    elseif s === :ha
        return ha(c)
    elseif s === :λ
        return λ(c)
    elseif s === :β
        return β(c)
    elseif s === :eclon
        return eclon(c)
    elseif s === :eclat
        return eclat(c)
    elseif s === :sgl
        return sgl(c)
    elseif s === :sgb
        return sgb(c)
    
    else
        throw(ArgumentError("Coordinate has no property $s"))
    end
end