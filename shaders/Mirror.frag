#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
 vec4 wo = v_position - vec4(u_cam_pos, 0);
 vec4 wi = wo - 2 * dot(wo, v_normal / length(v_normal)) * v_normal / length(v_normal);
 
  out_color = texture(u_texture_cubemap, wi.xyz / length(wi.xyz));

  out_color.a = 1;
}
