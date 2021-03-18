module Materials

using Images
using FileIO

using ..GfxBase

export Material, Texture
export ShadingModel, Flat, Normal
export PhysicalShadingModel, Lambertian, BlinnPhong

export get_diffuse

## Texture struct definition ##
mutable struct Texture
    image_data::Array{Vec3,2}
    repeat::Bool
end

""" Texture constructor - loads image data from a file"""
function Texture(image_fn::String, repeat::Bool)
    image_data = load(image_fn)
    Texture(image_data, repeat)
end

## Material struct definition ##
mutable struct Material
    mirror_coeff::Float32
    transparent_coeff::Float32
    diffuse_color::Union{Vec3,Nothing}
end

""" Get the diffuse Vec3 of a material; if the material is textured,
provide the uv coordinate on the object of that material. """
function get_diffuse(material::Material, uv::Union{Vec2,Nothing})
    material.diffuse_color
end

""" Look up a texture value given the uv coordinates """
function get_texture_value(texture::Texture, uv::Vec2)
    width, height = size(texture.image_data)
    texture.image_data[trunc(Int, (1 - uv[2]) * (height - 1)) + 1,
                       trunc(Int, uv[1] * (width - 1)) + 1]
end

end # module Materials
