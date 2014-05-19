#version 150

in  vec4 vPosition;
in  vec4 vColor;
out vec4 color;

uniform mat4 arcball;
uniform mat4 zoom;

void main() 
{   
	color = vColor;
	gl_Position = zoom * arcball * vPosition;
} 
