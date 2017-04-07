#version 3.7;

global_settings{ assumed_gamma 1.0 }

#include "library.inc"

// Camera
camera {
  angle 10
  location <40.0 ,40.0 ,-40.0>
  right x * image_width/image_height
  look_at <0.0 , 0.0, 0.0>
}

// Lighting
light_source { <1500,2500,-2500> color rgb<1, 1, 1> }
sky_sphere {
  pigment {
    gradient <0,1,0>
    color_map{
      [0 color rgb<1,1,1>]
      [0.4 color rgb<1,1,1>]
      [0.6 color rgb<1,1,1>]
      [1.0 color rgb<1,1,1>]
    }
    scale 2
  }
}

// Rasterized plane below objects
#declare constant_raster_scale = 0.5;
#declare constant_raster_half_line = 0.0185;

plane {
  <0,1,0>, 0
  texture {
    pigment { color rgbt<1, 1, 1, 0.7> * 1.1 }
    finish {ambient 0.45 diffuse 0.85}
  }
  texture {
    Raster(constant_raster_scale, constant_raster_half_line)
    rotate<0,0,0> 
  }
  texture {
    Raster(constant_raster_scale, constant_raster_half_line)
    rotate<0,90,0>
  }
  rotate<0,0,0>
  translate<0, -5, 0>
  no_shadow
}

// Scene settings
#default {
  finish {ambient 0.15 diffuse 0.85}
}
