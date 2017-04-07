#ifndef(SHORT_LIBRARY_INC)
#declare SHORT_LIBRARY_INC = version;
#version 3.6;

#ifdef(View_POV_Include_Stack)
     #debug "including analytical_g.inc\n"
#end

#include "math.inc"
#include "transforms.inc"
#include "functions.inc"

// Constant numeric declarations
#declare constant_atom_sphere_radius = 0.1;
#declare constant_bond_cylinder_radius = 0.05;
#declare constant_gradient_vector_parameter = 0.05;
#declare constant_tetrahedron_corner_sphere_radius = 0.15;
#declare constant_tetrahedron_edge_cylinder_radius = 0.1;

// Color declarations
#declare color_gray_ten = rgb<0.1, 0.1, 0.1>;
#declare color_gray_sixty = rgb<0.6, 0.6, 0.6>;
#declare color_steel_blue = rgb<0.274, 0.51, 0.7>;
#declare color_tomato = rgb<1.0, 0.38, 0.28>;

// Macros
#macro Vector(point_start,point_end, vector_radius)
  union{
    cylinder {
      point_start,
      point_end - ( vnormalize(point_end - point_start) * 9.5 * vector_radius),
      vector_radius
    }
    cone {
      point_end - ( vnormalize(point_end - point_start) * 10 * vector_radius),
      3 * vector_radius,
      point_end,
      0
    }
  }
#end

#macro Tetrahedron(corner_one, corner_two, corner_three, corner_four, sphere_radius, edge_cylinder_radius, fill)
  union {
    #if (sphere_radius != 0)
      sphere { corner_one, sphere_radius}
      sphere { corner_two, sphere_radius}
      sphere { corner_three, sphere_radius}
      sphere { corner_four, sphere_radius}
    #end

    #if (edge_cylinder_radius != 0)
      cylinder { corner_one, corner_two, edge_cylinder_radius }
      cylinder { corner_one, corner_three, edge_cylinder_radius }
      cylinder { corner_one, corner_four, edge_cylinder_radius }
      cylinder { corner_two, corner_three, edge_cylinder_radius }
      cylinder { corner_three, corner_four, edge_cylinder_radius }
      cylinder { corner_four, corner_two, edge_cylinder_radius }
    #end

    #if (fill = 1)
      triangle { corner_one, corner_two, corner_three }
      triangle { corner_one, corner_three, corner_four }
      triangle { corner_one, corner_four, corner_two }
      triangle { corner_four, corner_two, corner_three }
    #end
  }
#end

#macro Raster(scale_factor, half_line) 
  pigment { 
    gradient x scale scale_factor
    color_map { 
      [0.000 color rgbt<1,1,1,1> * 0.2]
      [0 + half_line color rgbt<1,1,1,1> * 0.2]
      [0 + half_line color rgbt<1,1,1,1>]
      [1 - half_line color rgbt<1,1,1,1>]
      [1 - half_line color rgbt<1,1,1,1> * 0.2]
      [1.000 color rgbt<1,1,1,1> * 0.2]
    }
  }
  finish { ambient 0.15 diffuse 0.85}
#end


// Composite Macros, for representing a Molecule
#macro Atom(point)
  object {
    sphere { point, constant_atom_sphere_radius}
    pigment { color color_gray_sixty }
  }
#end

#macro Bond(point_one, point_two)
  object {
    cylinder { point_one, point_two, constant_bond_cylinder_radius }
    pigment { color color_gray_sixty }
  }
#end

#macro GradientVector(point, gradientVector)
  #if (vlength(gradientVector) < 0.6)
    object {
      cylinder { point, point + gradientVector, constant_gradient_vector_parameter }
      pigment { color color_tomato }
    }
  #else
    object {
      Vector(point, point + gradientVector, constant_gradient_vector_parameter)
      pigment { color color_tomato }
    }
  #end
#end

#macro TetrahedronHighlight(corner_one, corner_two, corner_three, corner_four)
  object{
    Tetrahedron(
      corner_one,
      corner_two,
      corner_three,
      corner_four,
      constant_tetrahedron_corner_sphere_radius,
      constant_tetrahedron_edge_cylinder_radius,
      0
    )
    pigment { color color_steel_blue }
  }
#end

#version SHORT_LIBRARY_INC;
#end // ifndef
