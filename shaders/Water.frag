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


  vec3 h2 = normalize(u_cam_pos - v_position.xyz) + normalize(u_light_pos - v_position.xyz);

  h2 = h2 / length(h2);

  vec3 pos = (u_light_pos - v_position.xyz);

  vec3 ka = vec3(0.3, 0.3, 0.3);
  vec3 ks = vec3(0.5, 0.5, 0.5);
  vec3 kd = u_color.xyz;

//0, 0, 1 => Blue
//0, 1, 0 => Green
//1, 0, 0 => Red
  vec3 Ia = vec3(0, 0.8, 0);

  vec3 ambient = ka * Ia;

  vec3 I = u_light_intensity / dot(pos, pos);

  float p = 50.0;

  vec3 specular = ks * I * pow(clamp(dot(normalize(v_normal.xyz), h2), 0, 1), p);
  
  //this is color
  vec3 diffuse = kd * I * max(0, dot(normalize(v_normal.xyz), normalize(pos))) * u_color.xyz;

   out_color = vec4(diffuse + ambient + specular, 1);
  
  out_color.a = 1;
}

