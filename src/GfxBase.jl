module GfxBase

export Vec4, Vec3, Vec2, Ray, distance

using Images
using StaticArrays
import LinearAlgebra.cross
import LinearAlgebra.normalize

export RGB

# some type aliases:
const Vec4 = SVector{4, Float32} # 4-vector of floats
const Vec3 = SVector{3, Float32} # 3-vector of floats
const Vec2 = SVector{2, Float32} # 2-vector of floats

struct Ray
    origin::Vec3
    direction::Vec3
end

function distance(a::Vec3, b::Vec3)
    sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2)
end

end # module GfxBase
