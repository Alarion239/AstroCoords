module AstroCoords


include("transformations.jl")

export Coordinate, AbstractFrame, AbstractRepresentation, AbstractCartesianRepresentation, AbstractSphericalRepresentation, Cartesian, CartesianD, Spherical, SphericalD, ICRS, FK4, FK5, FK4NoETerms, HADec, AltAz, Ecliptic, Galactic, Supergalactic, TEME, CIRS, ITRS, HCRS

end
