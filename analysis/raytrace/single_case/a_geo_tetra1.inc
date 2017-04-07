// POV-Ray 3.7 Scene File "a_geo_tetra1.pov"
// Author: Friedrich A. Lohmueller, 2005/2009/2011

// email: Friedrich.Lohmueller_at_t-online.de
// homepage: http://wwww.f-lohmueller.de
//
//------------------------------------------------------------------------
#version 3.7;
global_settings{ assumed_gamma 1.0 }
//------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"          
#include "transforms.inc"

//------------------------------------------------------------------------
#declare Camera_0 = camera {                                 // xy-view
                             angle 20
                             location  <0.0 , 0.0 ,-38.0>
                             right     x*image_width/image_height
                             look_at   <0.0, 0.0 , 0.0>
                           }

#declare Camera_1 = camera {                                // diagonal view
                             angle 8
                             location  <40.0 ,40.0 ,-40.0>
                             right     x*image_width/image_height
                             look_at   <0.0 , 0.0, 0.0>
                           }
#declare Camera_2 = camera {                               // yz-view
                             angle 20          
                             location  <31.0 , 7.0 ,1.0>
                             right     x*image_width/image_height
                             look_at   <1.5 , 1.4 , 1.0>
                           }
#declare Camera_3 = camera { 
                             angle 20
                             location  < 0.0 , 37.5 ,-0.0001>
                             right     x*image_width/image_height
                             look_at   <0 , 0 , 2.0>
                           }

camera{Camera_1}

//------------------------------------------------------------------------
// sun -------------------------------------------------------------------
light_source{<1500,2500,-2500> color White}
// sky -------------------------------------------------------------------
sky_sphere{ pigment{ gradient <0,1,0>
                     color_map{ [0   color rgb<1,1,1>         ]//White
                                [0.4 color rgb<1,1,1>]//~Navy
                                [0.6 color rgb<1,1,1>]//<0.14,0.14,0.56>]//~Navy
                                [1.0 color rgb<1,1,1>         ]//White
                              }
                     scale 2 }
           } // end of sky_sphere 
//------------------------------------------------------------------------

//------------------------------ the Axes --------------------------------
//------------------------------------------------------------------------
#macro Axis_( AxisLen, Dark_Texture,Light_Texture) 
 union{
    cylinder { <0,-AxisLen,0>,<0,AxisLen,0>,0.05
               texture{checker texture{Dark_Texture } 
                               texture{Light_Texture}
                       translate<0.1,0,0.1>}
             }
    cone{<0,AxisLen,0>,0.2,<0,AxisLen+0.7,0>,0
          texture{Dark_Texture}
         }
     } // end of union                   
#end // of macro "Axis()"
//------------------------------------------------------------------------
#macro AxisXYZ( AxisLenX, AxisLenY, AxisLenZ, Tex_Dark, Tex_Light)
//--------------------- drawing of 3 Axes --------------------------------
union{
#if (AxisLenX != 0)
 object { Axis_(AxisLenX, Tex_Dark, Tex_Light)   rotate< 0,0,-90>}// x-Axis
 text   { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf",  "x",  0.15,  0  texture{Tex_Dark} 
          scale 0.5 translate <AxisLenX+0.05,0.4,-0.10>}
#end // of #if 
#if (AxisLenY != 0)
 object { Axis_(AxisLenY, Tex_Dark, Tex_Light)   rotate< 0,0,  0>}// y-Axis
 text   { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf",  "y",  0.15,  0  texture{Tex_Dark}    
           scale 0.5 translate <-0.55,AxisLenY+0.20,-0.10>}
#end // of #if 
#if (AxisLenZ != 0)
 object { Axis_(AxisLenZ, Tex_Dark, Tex_Light)   rotate<90,0,  0>}// z-Axis
 text   { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf",  "z",  0.15,  0  texture{Tex_Dark}
               scale 0.5 translate <-0.55,0.1,AxisLenZ+0.10>}
#end // of #if 
} // end of union
#end// of macro "AxisXYZ( ... )"

//------------------------------------------------------------------------
#declare Texture_A_Dark  = texture {
                               pigment{color rgb<0.4,0.4,0.4>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }
#declare Texture_A_Light = texture { 
                               pigment{color rgb<1,1,1>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }

// object{ AxisXYZ( 6.0, 5.0, 6.0, Texture_A_Dark, Texture_A_Light)} 
//-------------------------------------------------- end of coordinate axes


// ground -----------------------------------------------------------------
//---------------------------------<<< settings of squered plane dimensions
#declare RasterScale = 1.0/2;
#declare RasterHalfLine  = 0.035/2;
#declare RasterHalfLineZ = 0.035/2;
//-------------------------------------------------------------------------
#macro Raster(RScale, HLine) 
       pigment{ gradient x scale RScale
                color_map{[0.000   color rgbt<1,1,1,1>*0.2]
                          [0+HLine color rgbt<1,1,1,1>*0.2]
                          [0+HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,1>*0.2]
                          [1.000   color rgbt<1,1,1,1>*0.2]} }
       finish { ambient 0.15 diffuse 0.85}
 #end// of Raster(RScale, HLine)-macro
//-------------------------------------------------------------------------


plane { <0,1,0>, 0    // plane with layered textures
        texture { pigment{color rgbt<1,1,1,0.7>*1.1}
                  finish {ambient 0.45 diffuse 0.85}}
        texture { Raster(RasterScale,RasterHalfLine ) rotate<0,0,0> }
        texture { Raster(RasterScale,RasterHalfLineZ) rotate<0,90,0>}
        rotate<0,0,0>
        translate<0, -5, 0>
        no_shadow
      }
//------------------------------------------------ end of squered plane XZ

// Scene settings
#default{ finish {ambient 0.15 diffuse 0.85} } //

// Includes
#include "analytical_g.inc"

// Some scene settings
#declare Rl  = 0.1; // radius of lines  - Radius der Strecken 
#declare Rp  = 0.15; // radius of points - Radius der Punkte 

#declare myGray = rgb<0.1, 0.1, 0.1>;
#declare lightGray = rgb<0.6, 0.6, 0.6>;
#declare steelBlue = rgb<0.274, 0.51, 0.7>;
#declare tomato = rgb<1.0, 0.38, 0.28>;

/*fog{ 
  fog_type 2
  distance  25
  color     White
  fog_offset -5
  fog_alt 5
  up <1, 1, -1>
}*/

/* object{ Show_Yxz( P, Rl/2) pigment{ color lightGray }}
object{ Show_Yxz( Q, Rl/2) pigment{ color lightGray }}
object{ Show_Yxz( R, Rl/2) pigment{ color lightGray }}
object{ Show_Yxz( S, Rl/2) pigment{ color lightGray }}*/

// Point reference (mostly unneccessary)
/*text { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", "P",0.1,0
       scale 0.5 rotate<20,-45,0> translate P+< 0.5, 0.3,-0.1> pigment{ color myGray } no_shadow }
text { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", "Q",0.1,0
       scale 0.5 rotate<20,-45,0> translate Q+< 0.0,0.3, 0.0>  pigment{ color myGray } no_shadow }
text { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", "R",0.1,0
       scale 0.5 rotate<20,-45,0> translate R+<-0.5,0.1, 0.0>  pigment{ color myGray } no_shadow }
text { ttf "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", "S",0.1,0
       scale 0.5 rotate<20,-45,0> translate S+< 0.0,0.3, 0.0>  pigment{ color myGray } no_shadow }*/

