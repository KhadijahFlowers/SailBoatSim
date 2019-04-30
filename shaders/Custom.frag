#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_3;
uniform vec2 u_texture_3_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;


void main() {
  //BUMP
  vec4 newNorm = normalize(v_normal);
  vec4 newTan = normalize(v_tangent);

  vec3 bitangent = cross(newNorm.xyz, newTan.xyz);

  mat3 TBN = mat3(newTan, normalize(bitangent), newNorm);

  vec2 w = v_uv + vec2(1.0 / u_texture_3_size[0], 0.0);
  vec2 h = v_uv + vec2(0.0, 1.0 / u_texture_3_size[1]);


  float dU = (u_normal_scaling * u_height_scaling * (texture(u_texture_3, w)[0])) - 
  (u_normal_scaling * u_height_scaling * (texture(u_texture_3, v_uv)[0]));
 
  float dV = (u_normal_scaling * u_height_scaling * (texture(u_texture_3, h)[0])) - 
  (u_normal_scaling * u_height_scaling * (texture(u_texture_3, v_uv)[0]));
  
  vec3 localNorm = vec3(-dU, -dV, 1.0);


  vec3 displace = TBN * localNorm;

  out_color = vec4(displace, 0);
  
  out_color.a = 1;

}

