### Main Program: parse args and call main() ###	
using ArgParse	

using WWURay	

""" parse_cmdline	
Parse the command line args to specify scene, camera, and image size """	
function parse_cmdline()	
    s = ArgParseSettings()	
    @add_arg_table s begin	
        "--scene", "-s"	
            help="scene"	
            arg_type=String	
            default="teapot"	
        "--camera", "-c"	
            help="camera"	
            arg_type=Int	
            default=1	
        "width"	
            help="image width"	
            arg_type=Int
        "height"	
            help="image height"	
            arg_type=Int	
        "--samples"
            help="number of samples"
            arg_type=Int
        "--denoise"
            help="denoise resulting image"
            arg_type=Bool
    end	

    return parse_args(s)	
end	

a = parse_cmdline()	

WWURay.main(a["scene"], a["width"], a["height"], a["camera"], a["samples"], a["denoise"])