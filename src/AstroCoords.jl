module AstroCoords

# Include types first - TransformGraph needs these
include("epochs.jl")
include("frames.jl")
include("representation.jl")

# Now include transformations that use the types
include("transformations/graph.jl")

# Then include modules that use both types and transformations
include("coordinate.jl")
include("coordinate_accessors.jl")

function __init__()
    build_paths!()   # one-shot build
end

export 
  AbstractStandardEpoch, J2000, B1950, B1900,
  Coordinate, transform,
  Cartesian, Spherical, SphericalD,
  AbstractRepresentation,
  x_coord, y_coord, z_coord, lat, lon, dist,
  distance, norm,
  AbstractFrame, ICRS, FK4, FK5, FK4NoETerms, HADec, AltAz, Ecliptic, Galactic, Supergalactic, TEME, CIRS, ITRS, HCRS,
  ra, dec, gal, gab, alt, az, ha, β, λ, eclon, eclat, sgl, sgb,
  @transformation, build_paths!, transform
end
