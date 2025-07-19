"""
    Coordinates

Core coordinate type definitions and basic operations.

This module defines the fundamental `Coordinate{Frame, Representation}` type that combines
a reference frame with a coordinate representation. It provides the foundation for all
coordinate operations in AstroCoords.jl.
"""

using ..AstroCoords: AbstractFrame, AbstractRepresentation, Cartesian

"""
    Coordinate{F<:AbstractFrame, R<:AbstractRepresentation}

A coordinate in a specific reference frame with a specific representation.

The `Coordinate` type is the fundamental building block of AstroCoords.jl. It combines
a reference frame (like ICRS, Galactic, etc.) with a coordinate representation 
(like spherical or Cartesian coordinates).

# Type Parameters
- `F`: The reference frame type (subtype of `AbstractFrame`)
- `R`: The coordinate representation type (subtype of `AbstractRepresentation`)

# Fields
- `frame::F`: The reference frame instance
- `representation::R`: The coordinate representation instance

# Examples
```julia-repl
julia> using AstroCoords

julia> # Create an ICRS coordinate with spherical representation
julia> coord = Coordinate(ICRS(), SphericalD(deg2rad(45), deg2rad(30), 100.0))

julia> # Access the frame and representation
julia> coord.frame
ICRS()

julia> coord.representation
SphericalD{Float64, Float64}(0.5235987755982988, 0.7853981633974483, 100.0)
```

# Performance Notes
The coordinate type is designed for maximum performance:
- Type parameters enable compile-time optimizations
- Field access is zero-cost
- Frame-specific accessors compile to direct field access
"""
struct Coordinate{F<:AbstractFrame, R<:AbstractRepresentation}
    frame::F
    representation::R
    function Coordinate(frame::F, representation::R) where {F<:AbstractFrame, R<:AbstractRepresentation}
        new{F, R}(frame, representation)
    end
end

"""
    transform(c::Coordinate, target::AbstractFrame) -> Coordinate

Transform a coordinate from its current frame to a target frame.

This high-level interface handles representation conversions automatically:
1. Converts any representation to Cartesian
2. Calls the low-level transform
3. Wraps the result back into a Coordinate

# Examples
```julia-repl
julia> using AstroCoords

julia> # Create a coordinate in ICRS
julia> coord = Coordinate(ICRS(), SphericalD(deg2rad(45), deg2rad(30), 100.0))

julia> # Transform to Galactic frame
julia> gal_coord = transform(coord, Galactic())

julia> gal_coord.frame
Galactic()
```
"""
function transform(c::Coordinate{F,R}, target::T) where {
    F<:AbstractFrame, R<:AbstractRepresentation, T<:AbstractFrame}
    
    # 1. Convert any representation into Cartesian
    cart = R <: Cartesian ? c.representation : Cartesian(c.representation)
    
    # 2. Call the 3-arg graph engine — gives back a Cartesian
    new_cart = transform(cart, typeof(c.frame), typeof(target))
    
    # 3. Wrap result into a Coordinate
    return Coordinate(target, new_cart)
end