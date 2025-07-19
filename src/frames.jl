abstract type AbstractFrame end
Base.Broadcast.broadcastable(x::AbstractFrame) = Ref(x)

# ================= Fundamental Frames =================
struct ICRS <: AbstractFrame end
struct Galactic <: AbstractFrame end
struct Supergalactic <: AbstractFrame end

# ================= Equinox-Dependent Frames =================
struct FK5{E} <: AbstractFrame end


struct FK4{E} <: AbstractFrame end

struct FK4NoETerms{E} <: AbstractFrame end


# ================= Ecliptic Frames =================

abstract type EclipticFrameOrigin end 
struct Geocentric <: EclipticFrameOrigin end
struct Heliocentric <: EclipticFrameOrigin end
struct Barycentric <: EclipticFrameOrigin end

abstract type EclipticEquinoxType end
struct MeanEquinox <: EclipticEquinoxType end
struct TrueEquinox <: EclipticEquinoxType end
struct IAU76Equinox <: EclipticEquinoxType end

struct Ecliptic{Origin <: EclipticFrameOrigin, Equinox <: EclipticEquinoxType, E, O} <: AbstractFrame
    obstime :: O           
    equinox :: E           
    function Ecliptic(obstime::O, equinox::E, origin::Origin, equinox_type::Equinox) where {O, E, Origin <: EclipticFrameOrigin, Equinox <: EclipticEquinoxType}
        new{Origin, Equinox, typeof(obstime), typeof(equinox)}(obstime, equinox, origin, equinox_type)
    end
end

const GeocentricMeanEcliptic = Ecliptic{Geocentric, MeanEquinox}
const GeocentricTrueEcliptic = Ecliptic{Geocentric, TrueEquinox}
const BarycentricMeanEcliptic = Ecliptic{Barycentric, MeanEquinox}
const BarycentricTrueEcliptic = Ecliptic{Barycentric, TrueEquinox}
const HeliocentricMeanEcliptic = Ecliptic{Heliocentric, MeanEquinox}
const HeliocentricTrueEcliptic = Ecliptic{Heliocentric, TrueEquinox}
const HeliocentricIAU76Ecliptic = Ecliptic{Heliocentric, IAU76Equinox}

# ================= Observer-dependent Frames =================

struct AltAz{E, L, T} <: AbstractFrame 
    obstime :: E           
    location::L         
    pressure::T     
    temperature::T   
    relative_humidity::T  
    obswl::T          
    function AltAz(obstime::E, location::L, pressure::T, temperature::T, relative_humidity::T, obswl::T) where {E, L, T}
        new{typeof(obstime), typeof(location), T}(obstime, location, pressure, temperature, relative_humidity, obswl)
    end
end

struct HADec{E, L} <: AbstractFrame 
    obstime :: E           
    location:: L         
    function HADec(obstime::E, location::L) where {E, L}
        new{typeof(obstime), typeof(location)}(obstime, location)
    end
end

# ================= Intermediate Reference Frames =================

struct GCRS{E, P, V} <: AbstractFrame 
    obstime :: E           
    obspos  :: P
    obsvel :: V
    function GCRS(obstime::E, obspos::P, obsvel::V) where {E, P, V}
        new{typeof(obstime), typeof(obspos), typeof(obsvel)}(obstime, obspos, obsvel)
    end
end

struct CIRS{E} <: AbstractFrame 
    obstime :: E           
    function CIRS(obstime::E) where E
        new{typeof(obstime)}(obstime)
    end
end

struct ITRS{E} <: AbstractFrame 
    obstime :: E           
    function ITRS(obstime::E) where E
        new{typeof(obstime)}(obstime)
    end
end

struct HCRS{E} <: AbstractFrame 
    obstime :: E           
    function HCRS(obstime::E) where E
        new{typeof(obstime)}(obstime)
    end
end

struct TEME{E} <: AbstractFrame 
    obstime :: E           
    function TEME(obstime::E) where E
        new{typeof(obstime)}(obstime)
    end
end

struct TETE{E} <: AbstractFrame 
    obstime :: E           
    function TETE(obstime::E) where E
        new{typeof(obstime)}(obstime)
    end
end

struct PrecessedGeocentric{E, O, P, V} <: AbstractFrame 
    equinox :: E           
    obstime :: O           
    obspos::P
    obsvel::V 
    function PrecessedGeocentric(equinox::E, obstime::O, obspos::P, obsvel::P) where {E, O, P}
        new{typeof(equinox), typeof(obstime), typeof(obspos), typeof(obsvel)}(equinox, obstime, obspos, obsvel)
    end
end

# ================= Galactocentric Frame =================

struct Galactocentric{T} <: AbstractFrame 
    galcen_coord::T
    galcen_distance::T
    galcen_v_sun::T
    z_sun::T
    roll::T
    function Galactocentric(galcen_coord::T, galcen_distance::T, galcen_v_sun::T, z_sun::T, roll::T) where T
        new{T}(galcen_coord, galcen_distance, galcen_v_sun, z_sun, roll)
    end
end

# ================= Local Standard of Rest Frames =================

struct LSR{T} <: AbstractFrame 
    v_bary::T
    function LSR(v_bary::T) where T
        new{T}(v_bary)
    end
end

struct GalacticLSR{T} <: AbstractFrame 
    v_bary::T   
    function GalacticLSR(v_bary::T) where T
        new{T}(v_bary)
    end
end

struct LSRK{T} <: AbstractFrame 
    v_bary::T
    function LSRK(v_bary::T) where T
        new{T}(v_bary)
    end
end

struct LSRD{T} <: AbstractFrame 
    v_bary::T
    function LSRD(v_bary::T) where T
        new{T}(v_bary)
    end
end


