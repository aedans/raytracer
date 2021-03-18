""" 
Main module for CS480/580 A2 raytracer. Contains core raytracing algrithm,
while referencing several other modules that encapsulate supporting
functionality.
"""

module WWURay

export main, test

using FileIO
using Images
using StaticArrays
using LinearAlgebra
using BM3D
using Wavelets
using OpenCL

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Materials.jl")
include("WWUMeshes.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")

using .GfxBase
using .Materials
using .Scenes

import .Cameras
import .TestScenes

function toCL(vec::Vec3)
    Vec4(vec[1], vec[2], vec[3], 0)
end

function toCL(int::Int)
    convert(Int32, int)
end

function toCL(bool::Bool)
    if bool
        toCL(1)
    else
        toCL(0)
    end
end

struct CLRay
    origin::Vec4
    direction::Vec4
end

function toCL(ray::Ray)
    CLRay(toCL(ray.origin), toCL(ray.direction))
end

struct CLMaterial
    mirror_coeff::Float32
    transparent_coeff::Float32
    diffuse_color::Vec4
end

function toCL(material::Material)
    CLMaterial(material.mirror_coeff, material.transparent_coeff, toCL(material.diffuse_color))
end

struct CLLight
    color::Vec4
    point::Vec4
    radius::Float32
end

function toCL(light::Light)
    CLLight(toCL(light.color), toCL(light.point), light.radius)
end

struct CLTriangle
    a::Vec4
    b::Vec4
    c::Vec4
    anormal::Vec4
    bnormal::Vec4
    cnormal::Vec4
    auv::Vec2
    buv::Vec2
    cuv::Vec2
    material::Int32
end

function toCL(triangle::Triangle)
    normals = triangle.normals
    if normals === nothing
        normals = [Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)]
    end

    uvs = triangle.uvs
    if uvs === nothing
        uvs = [Vec2(0, 0), Vec2(0, 0), Vec2(0, 0)]
    end

    CLTriangle(
        toCL(triangle.positions[1]),
        toCL(triangle.positions[2]),
        toCL(triangle.positions[3]),
        toCL(normals[1]),
        toCL(normals[2]),
        toCL(normals[3]),
        uvs[1],
        uvs[2],
        uvs[3],
        toCL(triangle.material) - 1)
end

struct CLBounds
    min::Vec4
    max::Vec4
    leaf::Int32
    a::Int32
    b::Int32
end

function toCL(bounds::Bounds)
    CLBounds(toCL(bounds.min), toCL(bounds.max), toCL(bounds.leaf), toCL(bounds.a) - 1, toCL(bounds.b) - 1)
end

function test()
    main("bunny")
    main("cornell")
    main("teapot", 640, 360)
    # main("cubes")
    # main("glass")
end

function write(name, vecs, denoise)
    canvas = map(x -> RGB{Float32}(x[1], x[2], x[3]), vecs)

    if denoise
        canvas = bm3d(canvas, 10/255)
    end

    clamp01!(canvas)
    save(File(format"PNG", "out/" * name * ".png"), colorview(RGB, canvas))
end

function main(name, width=480, height=480, camera=1, samples=64, step_size=4, denoise=false)
    scene = TestScenes.get_scene(name)
    camera = TestScenes.get_camera(camera, height, width)

    source_file = open("src/ray.cl", "r")
    source = read(source_file, String)
    close(source_file)

    device, ctx, queue = cl.create_compute_context()
    program = cl.Program(ctx, source=source) |> cl.build!
    kernel = cl.Kernel(program, "ray_kernel")

    steps = samples / step_size
    size = height * width * step_size
    colors = zeros(Vec4, (height, width, step_size))
    canvas = zeros(Vec4, height, width)

    colors_buff = cl.Buffer(Vec4, ctx, :w, size)
    materials_buff = cl.Buffer(CLMaterial, ctx, (:r, :copy), length(scene.materials), hostbuf=map(x -> toCL(x), scene.materials))
    lights_buff = cl.Buffer(CLLight, ctx, (:r, :copy), length(scene.lights), hostbuf=map(x -> toCL(x), scene.lights))
    triangles_buff = cl.Buffer(CLTriangle, ctx, (:r, :copy), length(scene.triangles), hostbuf=map(x -> toCL(x), scene.triangles))
    bounds_buff = cl.Buffer(CLBounds, ctx, (:r, :copy), length(scene.bounds), hostbuf=map(x -> toCL(x), scene.bounds))
    
    @time for step in 1:steps
        @show step * step_size

        rays::Array{Ray, 1} = []
        for sample in 1:step_size, j in 1:width, i in 1:height
            push!(rays, Cameras.pixel_to_ray(camera, i + rand(Float64) - 0.5, j + rand(Float64) - 0.5))
        end

        seeds_buff = cl.Buffer(UInt32, ctx, (:r, :copy), size, hostbuf=rand(UInt32, size))
        rays_buff = cl.Buffer(CLRay, ctx, (:r, :copy), size, hostbuf=map(x -> toCL(x), rays))

        queue(kernel, size, nothing, 
            seeds_buff,
            rays_buff,
            materials_buff,
            lights_buff,
            toCL(length(scene.lights)),
            triangles_buff, 
            toCL(length(scene.triangles)),
            bounds_buff,
            toCL(length(scene.bounds)),
            scene.background,
            colors_buff)

        colors += reshape(cl.read(queue, colors_buff), (height, width, step_size))

        canvas = zeros(Vec4, height, width)
        for i in 1:height, j in 1:width, sample in 1:step_size
            canvas[i, j] += colors[i, j, sample] / (step_size * step)
        end

        write(name, canvas, false)
    end

    write(name, canvas, denoise)

    nothing
end

end 