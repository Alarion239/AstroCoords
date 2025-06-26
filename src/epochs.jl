abstract type AbstractStandardEpoch end

struct EpochJ2000 <: AbstractStandardEpoch end
struct EpochB1950 <: AbstractStandardEpoch end
struct EpochB1900 <: AbstractStandardEpoch end

const J2000 = EpochJ2000()
const B1950 = EpochB1950()
const B1900 = EpochB1900()