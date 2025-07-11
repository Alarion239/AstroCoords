module AstroCoords

# Include types first, then accessors that use those types
include("epochs.jl")
include("frames.jl")
include("representation.jl")
include("coordinate.jl")
include("coordinate_accessors.jl")

export 
  AbstractStandardEpoch, J2000, B1950, B1900,
  Coordinate,
  Cartesian, Spherical,
  AbstractRepresentation,
  x_coord, y_coord, z_coord, lat, lon, dist,
  distance, norm,
  AbstractFrame, ICRS, FK4, FK5, FK4NoETerms, HADec, AltAz, Ecliptic, Galactic, Supergalactic, TEME, CIRS, ITRS, HCRS,
  ra, dec, gal, gab, alt, az, ha, β, λ, eclon, eclat, sgl, sgb
end
