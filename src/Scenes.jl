module Scenes

export HitRecord, Scene, TriangleMesh, Triangle, Light, bvh, Bounds

using LinearAlgebra
using StaticArrays

using ..GfxBase
using ..WWUMeshes
using ..Materials

struct Triangle
    positions::SVector{3, Vec3}
    normals::Union{Nothing, SVector{3, Vec3}}
    uvs::Union{Nothing, SVector{3, Vec2}}
    material::Int
end

struct Bounds
    min::Vec3
    max::Vec3
    leaf::Bool
    a::Int
    b::Int
end

struct Light
    color::Vec3
    point::Vec3
    radius::Float32
end

"""Generic scene data type."""
struct Scene
    background::Vec3
    materials::Array{Material, 1}
    triangles::Array{Triangle, 1}
    lights::Array{Light, 1}
    bounds::Array{Bounds, 1}
end

function bounding_box(structures)  
    bb_min = [Inf, Inf, Inf]
    bb_max = [-Inf, -Inf, -Inf]

    for s in structures
        smin, smax = bounding_box(s)
        
        for i in 1:3
            bb_min[i] = min(bb_min[i], smin[i])
            bb_max[i] = max(bb_max[i], smax[i])
        end
    end

    Vec3(bb_min[1], bb_min[2], bb_min[3]), Vec3(bb_max[1], bb_max[2], bb_max[3])
end

function bounding_box(triangle::Triangle)
    t_min = [Inf, Inf, Inf]
    t_max = [-Inf, -Inf, -Inf]

    for n in 1:3
        p = triangle.positions[n]
        for i in 1:3
            t_min[i] = min(t_min[i], p[i])
            t_max[i] = max(t_max[i], p[i])
        end
    end

    Vec3(t_min[1], t_min[2], t_min[3]), Vec3(t_max[1], t_max[2], t_max[3])
end

function bounding_box(bounds::Bounds)
    bounds.min, bounds.max
end

function center_of(structure)
    min, max = bounding_box(structure)
    min .+ (max .- min) ./ 2
end

function bound_triangles(triangles, a::Int, b::Int)
    min, max = bounding_box([triangles[a], triangles[b]])
    Bounds(min, max, true, a, b)
end

function bound_bounds(bounds, a::Int, b::Int)
    min, max = bounding_box([bounds[a], bounds[b]])
    Bounds(min, max, false, a, b)
end

function bvh(triangles)
    bounds::Array{Bounds, 1} = []
    bvh(triangles, bounds, [i for i in 1:length(triangles)])
    bounds
end

function bvh(triangles, bounds::Array{Bounds, 1}, indices::Array{Int, 1})
    if length(indices) == 1
        push!(bounds, bound_triangles(triangles, indices[1], indices[1]))
        length(bounds)
    elseif length(indices) == 2
        push!(bounds, bound_triangles(triangles, indices[1], indices[2]))
        length(bounds)
    else
        i, j = farthest_pair(triangles, indices)

        is::Array{Int, 1} = []
        js::Array{Int, 1} = []

        for index in indices
            center = center_of(triangles[index])
            if distance(center_of(triangles[indices[i]]), center) < distance(center_of(triangles[indices[j]]), center)
                push!(is, index)
            else
                push!(js, index)
            end
        end

        ibounds = bvh(triangles, bounds, is)
        jbounds = bvh(triangles, bounds, js)
        push!(bounds, bound_bounds(bounds, ibounds, jbounds))
        length(bounds)
    end
end

function farthest_pair(triangles, indices::Array{Int, 1})
    max_dist = 0
    farthest_pair = nothing
    for i in 1:length(indices) - 1
        for j in i + 1:length(indices)
            dist = distance(center_of(triangles[indices[i]]), center_of(triangles[indices[j]]))
            if dist > max_dist
                max_dist = dist
                farthest_pair = i, j
            end
        end
    end

    farthest_pair
end

end