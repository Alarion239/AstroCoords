module UnitfulExt

using AstroCoords
using Unitful

# TODO: Implement Unitful support
Spherical(::T, ::T) where {T <: Unitful.Quantity} = throw(ArgumentError("UNITFUL IS NOT SUPPORTED"))
SphericalD(::T, ::T) where {T <: Unitful.Quantity} = throw(ArgumentError("UNITFUL IS NOT SUPPORTED"))
Cartesian(::T, ::T) where {T <: Unitful.Quantity} = throw(ArgumentError("UNITFUL IS NOT SUPPORTED"))
CartesianD(::T, ::T) where {T <: Unitful.Quantity} = throw(ArgumentError("UNITFUL IS NOT SUPPORTED"))  

# TODO: ensure that only the corresponding unitful type is used for the representation

end