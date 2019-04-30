#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;


void main() {
  // YOUR CODE HERE
  //vector in object space => transform back to model TBN Matrix
  // (Placeholder code. You will want to replace it.)

  //v_normal = normalize(u_model * in_normal);

  vec4 newNorm = normalize(v_normal);
  vec4 newTan = normalize(v_tangent);

  vec3 bitangent = cross(newNorm.xyz, newTan.xyz);

  mat3 TBN = mat3(newTan, normalize(bitangent), newNorm);

  vec2 w = v_uv + vec2(1.0 / u_texture_2_size[0], 0.0);
  vec2 h = v_uv + vec2(0.0, 1.0 / u_texture_2_size[1]);


  float dU = (u_normal_scaling * u_height_scaling * (texture(u_texture_2, w)[0])) - 
  (u_normal_scaling * u_height_scaling * (texture(u_texture_2, v_uv)[0]));
 
  float dV = (u_normal_scaling * u_height_scaling * (texture(u_texture_2, h)[0])) - 
  (u_normal_scaling * u_height_scaling * (texture(u_texture_2, v_uv)[0]));
  
  vec3 localNorm = vec3(-dU, -dV, 1.0);


  vec3 displace = TBN * localNorm;


  vec3 h2 = normalize(u_cam_pos - displace) + normalize(u_light_pos - displace);

  h2 = h2 / length(h2);

  vec3 pos = (u_light_pos - displace);

  vec3 ka = vec3(0.1, 0.1, 0.1);
  vec3 ks = vec3(0.5, 0.5, 0.5);
  vec3 kd = u_color.xyz;

  vec3 Ia = vec3(1, 1, 1);

  vec3 ambient = ka * Ia;

  vec3 I = u_light_intensity / dot(pos, pos);

  float p = 100.0;

  vec3 specular = ks * I * pow(clamp(dot(normalize(v_normal.xyz), h2), 0, 1), p);
  
  //this is color
  vec3 diffuse = kd * I * max(0, dot(normalize(v_normal.xyz), normalize(pos))) * u_color.xyz;

   out_color = vec4(ambient + diffuse + specular, 0);
  
  out_color.a = 1;
}

