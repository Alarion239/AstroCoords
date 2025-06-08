module AstroCoords

# Include types first, then accessors that use those types
include("frames.jl")
include("representation.jl")
include("coordinate.jl")
include("coordinate_accessors.jl")

export Coordinate,
  Cartesian, CartesianD, Spherical, SphericalD,
  AbstractRepresentation, AbstractCartesianRepresentation, AbstractSphericalRepresentation,
  x_coord, y_coord, z_coord, lat, lon, dist,
  distance, norm,
  AbstractFrame, ICRS, FK4, FK5, FK4NoETerms, HADec, AltAz, Ecliptic, Galactic, Supergalactic, TEME, CIRS, ITRS, HCRS

end
