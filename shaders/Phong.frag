#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
  //out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;


  vec3 h2 = u_light_pos + u_cam_pos;

 

  vec3 pos = (u_light_pos - v_position.xyz);

  float ka = 0.3;
  float ks = 0.5;
  float kd = 0.2;

  vec3 Ia = vec3(1, 1, 1);

  vec3 ambient = ka * Ia;

  vec3 I = u_light_intensity / (length(pos) * length(pos));

  float p = 100.0;

  vec3 specular = ks * I * pow(max(dot(v_normal.xyz, normalize(h2)), 0), p);

  vec3 diffuse = kd * I * max(0, dot(v_normal.xyz, normalize(pos)));

  vec3 res = diffuse + ambient + specular;
  // + ambient + specular;

  //diffuse;
  // + specular + ambient;
  
  out_color = vec4(res, 0);
  out_color.a = 1;
}

