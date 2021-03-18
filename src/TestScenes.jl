module TestScenes

using ..GfxBase
using ..Scenes
using ..Materials
using ..WWUMeshes
using ..Cameras

make_diffuse(color) = Material(0.0, 0.0, color)

black = Vec3(0,0,0)
red = Vec3(1,0,0)
green = Vec3(0,1,0)
blue = Vec3(0,0,1)
white = Vec3(1,1,1)
purple = Vec3(1,0,1)

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    triangles = []

    for t in mesh.triangles
        positions = [mesh.positions[t.positions[1]], mesh.positions[t.positions[2]], mesh.positions[t.positions[3]]]

        normals = nothing
        if t.normals !== nothing && length(t.normals) == 3
            normals = [mesh.normals[t.normals[1]], mesh.normals[t.normals[2]], mesh.normals[t.normals[3]]]
        end

        uvs = nothing
        if t.uvs !== nothing && length(t.uvs) == 3
            uvs = [mesh.uvs[t.uvs[1]], mesh.uvs[t.uvs[2]], mesh.uvs[t.uvs[3]]]
        end

        push!(triangles, Triangle(positions, normals, uvs, material))
    end

    triangles
end

function mesh_helper(mesh, material, scale=Vec3(1.0,1.0,1.0), translation=Vec3(0,0,0), rotation=Vec3(0, 0, 0))
    Rx = [1 0                0
          0 cos(rotation[1]) -sin(rotation[1])
          0 sin(rotation[1]) cos(rotation[1])]

    Ry = [cos(rotation[2])  0 sin(rotation[2])
          0                 1 0
          -sin(rotation[2]) 0 cos(rotation[2])]

    Rz = [cos(rotation[3]) -sin(rotation[3]) 0
          sin(rotation[3]) cos(rotation[3])  0
          0                0                 1]

    R = Rx * Ry * Rz

    for i in 1:length(mesh.positions)
        mesh.positions[i] = (R * mesh.positions[i]) .* scale + translation
    end

    create_triangles(mesh, material)
end

function camera_1(img_height, img_width)
    CanonicalCamera(img_height, img_width)
end

function camera_2(img_height, img_width)
    eye = Vec3(20, 4, 10)
    view = Vec3(-1, 0, -5) - eye
    up = Vec3(0, 1, 0)
    focal = 8.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end

function camera_3(img_height, img_width)
    Cameras.PerspectiveCamera(Vec3(-1, 0.8, -1.2), Vec3(1, -1, -1), Vec3(0, 1, 0), 0.3, img_height, img_width)
end

cameras = [camera_1, camera_2, camera_3]

function get_camera(i, img_height, img_width)
    cameras[i](img_height, img_width)
end

function bunny()
    bunny = read_obj("data/bunny.obj")
    bunny = mesh_helper(bunny, 1, 1.0, Vec3(0.2, 0, -5))

    cube = mesh_helper(cube_mesh(), 2, 10.0, Vec3(-11.2, 0, 0))

    lights = [ Light(Float32(1.3) * white, Vec3(1,2,-5), 1),
               Light(Float32(0.8) * white, Vec3(0,0,10), 3),
               Light(Float32(0.8) * white, Vec3(0,10,10), 3),
               Light(Float32(0.8) * white, Vec3(10,10,10), 3),
               Light(Float32(0.8) * white, Vec3(0,10,0), 3) ]

    materials = [Material(0.0, 0.0, Vec3(0.6, 0.5, 0.5)), Material(0.6, 0.0, white)]
    triangles = vcat(cube, bunny)
    Scene(black, materials, triangles, lights, bvh(triangles))
end

function teapot()
    teapot = read_obj("data/teapot.obj")
    estimate_normals(teapot)
    teapot = mesh_helper(teapot, 1, 0.5, Vec3(0, -1, -5))

    cube = mesh_helper(cube_mesh(), 1, 10.0, Vec3(0, -11, 0))

    lights = [ Light(4 * white, Vec3(1,2,-5), 1) ]

    materials = [Material(0.0, 0.0, white)]
    triangles = vcat(cube, teapot)
    Scene(black, materials, triangles, lights, bvh(triangles))
end

function cornell()
    cube_mat = Material(0, 0.0, .7 * white)

    cube_back = mesh_helper(cube_mesh(), 1, 2.0, Vec3(0, 0, -8))
    cube_left = mesh_helper(cube_mesh(), 2, 2.0, Vec3(-4, 0, -6))
    cube_right = mesh_helper(cube_mesh(), 3, 2.0, Vec3(4, 0, -6))
    cube_top = mesh_helper(cube_mesh(), 1, 2.0, Vec3(0, 4, -6))
    cube_bottom = mesh_helper(cube_mesh(), 1, 2.0, Vec3(0, -4, -6))

    cube_1 = mesh_helper(cube_mesh(), 1, Vec3(0.5, 1.0, 0.5), Vec3(-0.75, -1, -5), Vec3(0, pi / 8, 0))
    cube_2 = mesh_helper(cube_mesh(), 1, 0.5, Vec3(0.75, -1.5, -5), Vec3(0, -pi / 8, 0))

    lights = [ Light(4 * white, Vec3(0,1.49,-4), 0.5) ]

    materials = [Material(0, 0.0, .7 * white), Material(0, 0.0, .7 * red), Material(0, 0.0, .7 * green)]
    triangles = vcat(cube_back, cube_left, cube_right, cube_top, cube_bottom, cube_1, cube_2)
    Scene(.7 * white, materials, triangles, lights, bvh(triangles))
end

function cubes()
    cube_bottom = mesh_helper(cube_mesh(), make_diffuse(white), 6.0, Vec3(0, -9, -12))
    cube_back = mesh_helper(cube_mesh(), make_diffuse(blue), 6.0, Vec3(0, 0, -16))
    cube_opaque = mesh_helper(cube_mesh(), make_diffuse(green), 0.5, Vec3(-2, -1, -7))
    cube_semi_transparent = mesh_helper(cube_mesh(), Material(0.0, 0.75, white), 0.5, Vec3(0, -1, -7))

    lights = [ Light(4 * white, Vec3(0, 10, 10), 3)]
    
    triangles = vcat(cube_bottom, cube_back, cube_opaque, cube_semi_transparent)
    Scene(black, triangles, lights, bvh(triangles))
end

function glass()
    glass_mat = Material(0, 0.5, white)
    glass = read_obj("data/glass.obj")
    estimate_normals(glass)
    glass = mesh_helper(glass, glass_mat, 0.01, Vec3(0, 0.1, -3))

    cube_mat = Material(0.6, 0.0, white)
    cube = mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(0, -11, 0))

    lights = [ Light(4 * white, Vec3(1,2,-2), 1) ]

    triangles = vcat(cube, glass)
    Scene(black, triangles, lights, bvh(triangles))
end

function get_scene(name)
    scenes = Dict("bunny" => bunny(), "cornell" => cornell(), "teapot" => teapot())
    scenes[name]
end

end