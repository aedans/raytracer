module WWUMeshes

export read_obj, write_obj
export gen_mesh, est_normals
export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals
export OBJTriangle, OBJMesh

using FileIO
using LinearAlgebra

push!(LOAD_PATH, pwd())

include("GfxBase.jl")
using .GfxBase


""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int, 1} # vertex position indices
    uvs::Array{Int, 1} # vertex texture coordinate indices
    normals::Array{Int, 1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3, 1} # all vertex positions
    uvs::Array{Vec2, 1} # all texture coordinates
    normals::Array{Vec3, 1} # all vertex normals
    triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
end

""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float32, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float32, x) for x in tokens[2:3]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float32, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                if length(tokens) == 4
                    # create a OBJTriangle face:
                    points = []
                    uvs = []
                    normals = []
                    # handle faces with no texture and/or normals
                    for corner in tokens[2:4]
                        indices = split(corner, '/')
                        if length(indices) == 3 # all 3 present, third is a normal
                            push!(normals, parse(Int, indices[3]))
                        end
                        if length(indices) >= 2 && indices[2] != ""
                            # if there are 2 or more and the second isn't blank, it's a texture
                            push!(uvs, parse(Int, indices[2]))
                        end
                        if length(indices) >= 1 # first value is the position
                            push!(points, parse(Int, indices[1]))
                        else # unless it has none, in which case it's not valid
                            error("in line $line_number: face vertex $corner could not be parsed")
                        end
                    end
                    # create the triangle and add it to the triangles array
                    push!(m.triangles, OBJTriangle(points, uvs, normals))
                else
                    # create a OBJTriangle face:
                    points = []
                    uvs = []
                    normals = []
                    # handle faces with no texture and/or normals
                    for corner in tokens[2:4]
                        indices = split(corner, '/')
                        if length(indices) == 3 # all 3 present, third is a normal
                            push!(normals, parse(Int, indices[3]))
                        end
                        if length(indices) >= 2 && indices[2] != ""
                            # if there are 2 or more and the second isn't blank, it's a texture
                            push!(uvs, parse(Int, indices[2]))
                        end
                        if length(indices) >= 1 # first value is the position
                            push!(points, parse(Int, indices[1]))
                        else # unless it has none, in which case it's not valid
                            error("in line $line_number: face vertex $corner could not be parsed")
                        end
                    end
                    # create the triangle and add it to the triangles array
                    push!(m.triangles, OBJTriangle(points, uvs, normals))
                    
                    points = []
                    uvs = []
                    normals = []
                    # handle faces with no texture and/or normals
                    for corner in [tokens[4], tokens[5], tokens[2]]
                        indices = split(corner, '/')
                        if length(indices) == 3 # all 3 present, third is a normal
                            push!(normals, parse(Int, indices[3]))
                        end
                        if length(indices) >= 2 && indices[2] != ""
                            # if there are 2 or more and the second isn't blank, it's a texture
                            push!(uvs, parse(Int, indices[2]))
                        end
                        if length(indices) >= 1 # first value is the position
                            push!(points, parse(Int, indices[1]))
                        else # unless it has none, in which case it's not valid
                            error("in line $line_number: face vertex $corner could not be parsed")
                        end
                    end
                    # create the triangle and add it to the triangles array
                    push!(m.triangles, OBJTriangle(points, uvs, normals))
                end
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end

    end

end

""" tri_vertex_str(triangle)
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal) """
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end


""" gen_mesh(outfile, geom, divisionsU, divisionsV)
Generate a mesh and save the result in a file with name outfile.
geom may be "cube", "cylinder", or "sphere".
Cylinder requires divisionsU; sphere requires divisionsU and divisionsV. """
function gen_mesh(outfile, geom, divisionsU=0, divisionsV=0)
    if geom == "cube"
        mesh = cube_mesh()
    elseif geom == "cylinder"
        mesh = cylinder_mesh(divisionsU)
    elseif geom == "sphere"
        mesh = sphere_mesh(divisionsU, divisionsV)
    end
    write_obj(outfile, mesh)
end


""" est_normals(outfile, infile)
Estimate normals of the mesh stored in infile, saving the result in outfile."""
function est_normals(outfile, infile)
    input_mesh = read_obj(infile)
    mesh = estimate_normals(input_mesh)
    write_obj(outfile, mesh)
end


""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3( 1, -1, -1)) # 1 RBC
    push!(positions, Vec3( 1, -1,  1)) # 2 RBF
    push!(positions, Vec3(-1, -1,  1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3( 1,  1, -1)) # 5 RTC
    push!(positions, Vec3( 1,  1,  1)) # 6 RTF
    push!(positions, Vec3(-1,  1,  1)) # 7 LTF
    push!(positions, Vec3(-1,  1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3( 1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3( 0, 1, 0)) # U
    push!(normals, Vec3( 0,-1, 0)) # D
    push!(normals, Vec3( 0, 0, 1)) # C
    push!(normals, Vec3( 0, 0,-1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
    push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
    push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
    push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
    push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
    push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
    push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
    push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
    push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
    push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
    push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
    push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane. See the
assignment writeup for a diagram and details.
"""
function cylinder_mesh(divisionsU)
    positions = []
    uvs = []
    normals = []
    triangles = []

    # bottom center vertex
    push!(positions, Vec3(0, -1, 0))
    push!(uvs, Vec2(0.25, 0.75))
    push!(normals, Vec3(0, -1, 0))
    
    # top center vertex
    push!(positions, Vec3(0, 1, 0))
    push!(uvs, Vec2(0.75, 0.75))
    push!(normals, Vec3(0, 1, 0))

    # partially unwind loop where i = 0
    push!(positions, Vec3(1, 1, 0))
    push!(positions, Vec3(1, -1, 0))
    
    push!(uvs, Vec2(0.5, 0.75))
    push!(uvs, Vec2(1, 0.75))
    push!(uvs, Vec2(1, 0.5))
    push!(uvs, Vec2(1, 0))
    
    push!(normals, Vec3(1, 0, 0))
    
    delta = (2.0 * pi) / divisionsU
    for i in 1:divisionsU
        theta = delta * i
        x = cos(theta)
        z = sin(theta)

        push!(positions, Vec3(x, 1, z))
        push!(positions, Vec3(x, -1, z))

        push!(uvs, Vec2(0.25 + x / 4, 0.75 + z / 4))
        push!(uvs, Vec2(0.75 + x / 4, 0.75 + z / 4))
        push!(uvs, Vec2(1 - theta / (2.0 * pi), 0.5))
        push!(uvs, Vec2(1 - theta / (2.0 * pi), 0))

        push!(normals, Vec3(x, 0, z))

        # sides
        push!(triangles, OBJTriangle(
            [i * 2 + 2, i * 2 + 1, i * 2 + 3], 
            [i * 4 + 2, i * 4 + 1, i * 4 + 5],
            [i + 2, i + 2, i + 3]))

        push!(triangles, OBJTriangle(
            [i * 2 + 3, i * 2 + 4, i * 2 + 2], 
            [i * 4 + 5, i * 4 + 6, i * 4 + 2],
            [i + 3, i + 3, i + 2]))

        # caps
        push!(triangles, OBJTriangle(
            [1, i * 2 + 2, i * 2 + 4], 
            [1, i * 4 - 1, i * 4 + 3],
            [1, 1, 1]))

        push!(triangles, OBJTriangle(
            [i * 2 + 1, 2, i * 2 + 3], 
            [i * 4, 2, i * 4 + 4],
            [2, 2, 2]))
    end
    
    OBJMesh(positions, uvs, normals, triangles)
end


""" sphere_mesh(n, m)
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0), and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(n, m)
    positions = []
    uvs = []
    normals = []
    triangles = []

    deltaM = pi / m
    deltaN = (2.0 * pi) / n
    deltaV = 1 / m
    deltaU = 1 / n
    i = 0
    
    # partially unwind loop where j = 0
    for k in 1:n
        azimuth = deltaN * k

        push!(positions, Vec3(0, 1, 0))
        push!(uvs, Vec2(deltaU * k, 1))
        push!(normals, Vec3(0, 1, 0))
    end

    for j in 1:m + 1
        inclination = deltaM * j

        y = cos(inclination)
        
        # partially unwind loop where k = 0
        x = sin(inclination)

        push!(positions, Vec3(x, y, 0))
        push!(uvs, Vec2(0, 1 - deltaV * j))
        push!(normals, normalize(Vec3(x, y, 0)))

        for k in 1:n
            azimuth = deltaN * k

            x = sin(inclination) * cos(azimuth)
            z = sin(inclination) * sin(azimuth)

            i += 1
            push!(positions, Vec3(x, y, z))
            push!(uvs, Vec2(deltaU * k, 1 - deltaV * j))
            push!(normals, normalize(Vec3(x, y, z)))

            push!(triangles, OBJTriangle(
                [i, i + 1, i + n + 1], 
                [i, i + 1, i + n + 1], 
                [i, i + 1, i + n + 1]))
                
            push!(triangles, OBJTriangle(
                [i + n, i, i + n + 1], 
                [i + n, i, i + n + 1], 
                [i + n, i, i + n + 1]))
        end
    end

    OBJMesh(positions, uvs, normals, triangles)
end

"""
    estimate_normals(mesh::OBJMesh)
Estimates normals for the given mesh. Overwrites any existing normals and returns a new OBJMesh object.
"""
function estimate_normals(mesh::OBJMesh)
    mesh.normals = zeros(Vec3, length(mesh.positions))

    for t in 1:length(mesh.triangles)
        triangle = mesh.triangles[t]
        triangle.normals = zeros(3)

        v = mesh.positions[triangle.positions[2]] - mesh.positions[triangle.positions[1]]
        w = mesh.positions[triangle.positions[3]] - mesh.positions[triangle.positions[1]]
        
        normal = Vec3(
            (v[2] * w[3]) - (v[3] * w[2]),
            (v[3] * w[1]) - (v[1] * w[3]),
            (v[1] * w[2]) - (v[2] * w[1]))

        for v in 1:3
            mesh.normals[triangle.positions[v]] += normal
            triangle.normals[v] = triangle.positions[v]
        end
    end

    for i in 1:length(mesh.normals)
        mesh.normals[i] = normalize(mesh.normals[i])
    end

    mesh
end

end # module WWUMeshes


