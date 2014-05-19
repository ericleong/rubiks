// ECE462 Rubiks Cube
// Eric Leong
// 4/2/2013

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "Angel.h"

const int SCREEN_HEIGHT = 512;
const int SCREEN_WIDTH = 512;
const float TWIST_MIN_DIST = .08;
const float TWIST_MAX_ANGLE = 40; // degrees
const float TWIST_RATE = M_PI/180.0f/2.0f; // radians/pixel
const float SCALE_FACTOR = 2;

const float CUBELET_SIZE = 1;
const float CUBE_OFFSET = 2.15/2;
const float CUBE_SIZE = 2 * (CUBE_OFFSET + CUBELET_SIZE/2);

typedef Angel::vec4 color4;
typedef Angel::vec4 point4;

const int NumVertices = 36 * 27; //(6 faces)(2 triangles/face)(3 vertices/triangle)(27 cubes)

point4 oPoints[NumVertices];
point4 points[NumVertices];
point4 oColors[NumVertices];
color4 colors[NumVertices];

const int Xaxis = 0;
const int Yaxis = 1;
const int Zaxis = 2;

// Vertices of a unit cube centered at origin, sides aligned with axes
point4 vertices[8] = {
		point4(-CUBELET_SIZE/2, -CUBELET_SIZE/2, CUBELET_SIZE/2, 1.0),
		point4(-CUBELET_SIZE/2, CUBELET_SIZE/2, CUBELET_SIZE/2, 1.0),
		point4(CUBELET_SIZE/2, CUBELET_SIZE/2, CUBELET_SIZE/2, 1.0),
		point4(CUBELET_SIZE/2, -CUBELET_SIZE/2, CUBELET_SIZE/2, 1.0),
		point4(-CUBELET_SIZE/2, -CUBELET_SIZE/2, -CUBELET_SIZE/2, 1.0),
		point4(-CUBELET_SIZE/2, CUBELET_SIZE/2, -CUBELET_SIZE/2, 1.0),
		point4(CUBELET_SIZE/2, CUBELET_SIZE/2, -CUBELET_SIZE/2, 1.0),
		point4(CUBELET_SIZE/2, -CUBELET_SIZE/2, -CUBELET_SIZE/2, 1.0)
};

// RGBA colors
color4 vertex_colors[7] = {
		color4(0.2, 0.2, 0.2, 1.0),  // dark gray
		color4(1.0, 1.0, 1.0, 1.0),  // white
		color4(1.0, 0.0, 0.0, 1.0),  // red
		color4(0.0, 0.0, 1.0, 1.0),  // blue
		color4(1.0, 1.0, 0.0, 1.0),  // yellow
		color4(1.0, 0.5, 0.0, 1.0),  // orange
		color4(0.0, 1.0, 0.0, 1.0),  // green
	};

GLfloat Arcball[16];
GLfloat Zoom[16];

GLuint arcball; // The location of the "arcball" shader uniform variable
GLuint zoom; // The location of the "zoom" shader uniform variable

float scale = 1.0f/3.0f;

// Mouse interaction
int oldX, oldY;
bool mouse1 = false;
bool mouse2 = false;

bool selected = false;
struct plane {
	Angel::vec3 point;
	Angel::vec3 normal;
	Angel::mat4 rotate;
	int rotation;
} selectedPlane;
Angel::vec3 selectStart;
struct twister {
	Angel::vec3 axis;
	GLfloat angle;
	int layer;
	int dir;
} twist;

Angel::mat4 mArcball = RotateX(45) * RotateY(45) * RotateZ(45);
Angel::mat4 mTwist = Angel::mat4();
Angel::mat4 mZoom =  Angel::Scale(scale, scale, scale);

// save
char* saveName;
char* loadName;
bool loadAFile;

// randomize
int turns = 0;

//----------------------------------------------------------------------------
// MATH FUNCTIONS

int sign(float x) {
	return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

bool isNumber(float x) {
	// This looks like it should always be true,
	// but it's false if x is a NaN.
	return (x == x);
}

Angel::vec3 orthoShift(Angel::vec3 orthoAxis) {
	return Angel::vec3(orthoAxis.z, orthoAxis.x, orthoAxis.y);
}

Angel::vec3 transform3(Angel::mat4 transform, Angel::vec3 v) {
	// apply transformation
	Angel::vec4 u = transform * Angel::vec4(v.x, v.y, v.z, 0);
	return Angel::vec3(u.x, u.y, u.z);
}

// analytic determinate of 4x4 matrix
// https://bitbucket.org/t3hprogrammer/geoar/src/8ddc3c35f08d8ad77b6947b4bff1a6002d15cf6d/GeoAR%20App/jni/SampleMath.cpp
float det4(Angel::mat4 m) {
	return m[0].w * m[1].z * m[2].y * m[3].x - m[0].z * m[1].w * m[2].y * m[3].x
			- m[0].w * m[1].y * m[2].z * m[3].x
			+ m[0].y * m[1].w * m[2].z * m[3].x
			+ m[0].z * m[1].y * m[2].w * m[3].x
			- m[0].y * m[1].z * m[2].w * m[3].x
			- m[0].w * m[1].z * m[2].x * m[3].y
			+ m[0].z * m[1].w * m[2].x * m[3].y
			+ m[0].w * m[1].x * m[2].z * m[3].y
			- m[0].x * m[1].w * m[2].z * m[3].y
			- m[0].z * m[1].x * m[2].w * m[3].y
			+ m[0].x * m[1].z * m[2].w * m[3].y
			+ m[0].w * m[1].y * m[2].x * m[3].z
			- m[0].y * m[1].w * m[2].x * m[3].z
			- m[0].w * m[1].x * m[2].y * m[3].z
			+ m[0].x * m[1].w * m[2].y * m[3].z
			+ m[0].y * m[1].x * m[2].w * m[3].z
			- m[0].x * m[1].y * m[2].w * m[3].z
			- m[0].z * m[1].y * m[2].x * m[3].w
			+ m[0].y * m[1].z * m[2].x * m[3].w
			+ m[0].z * m[1].x * m[2].y * m[3].w
			- m[0].x * m[1].z * m[2].y * m[3].w
			- m[0].y * m[1].x * m[2].z * m[3].w
			+ m[0].x * m[1].y * m[2].z * m[3].w;
}

// analytic inverse of 4x4 matrix
// https://bitbucket.org/t3hprogrammer/geoar/src/8ddc3c35f08d8ad77b6947b4bff1a6002d15cf6d/GeoAR%20App/jni/SampleMath.cpp
Angel::mat4 inverse4(Angel::mat4 m) {
	Angel::mat4 r;

	float det = 1.0f / det4(m);

	r[0].x = m[2].y * m[3].z * m[1].w - m[3].y * m[2].z * m[1].w
			+ m[3].y * m[1].z * m[2].w - m[1].y * m[3].z * m[2].w
			- m[2].y * m[1].z * m[3].w + m[1].y * m[2].z * m[3].w;

	r[0].y = m[3].x * m[2].z * m[1].w - m[2].x * m[3].z * m[1].w
			- m[3].x * m[1].z * m[2].w + m[1].x * m[3].z * m[2].w
			+ m[2].x * m[1].z * m[3].w - m[1].x * m[2].z * m[3].w;

	r[0].z = m[2].x * m[3].y * m[1].w - m[3].x * m[2].y * m[1].w
			+ m[3].x * m[1].y * m[2].w - m[1].x * m[3].y * m[2].w
			- m[2].x * m[1].y * m[3].w + m[1].x * m[2].y * m[3].w;

	r[0].w = m[3].x * m[2].y * m[1].z - m[2].x * m[3].y * m[1].z
			- m[3].x * m[1].y * m[2].z + m[1].x * m[3].y * m[2].z
			+ m[2].x * m[1].y * m[3].z - m[1].x * m[2].y * m[3].z;

	r[1].x = m[3].y * m[2].z * m[0].w - m[2].y * m[3].z * m[0].w
			- m[3].y * m[0].z * m[2].w + m[0].y * m[3].z * m[2].w
			+ m[2].y * m[0].z * m[3].w - m[0].y * m[2].z * m[3].w;

	r[1].y = m[2].x * m[3].z * m[0].w - m[3].x * m[2].z * m[0].w
			+ m[3].x * m[0].z * m[2].w - m[0].x * m[3].z * m[2].w
			- m[2].x * m[0].z * m[3].w + m[0].x * m[2].z * m[3].w;

	r[1].z = m[3].x * m[2].y * m[0].w - m[2].x * m[3].y * m[0].w
			- m[3].x * m[0].y * m[2].w + m[0].x * m[3].y * m[2].w
			+ m[2].x * m[0].y * m[3].w - m[0].x * m[2].y * m[3].w;

	r[1].w = m[2].x * m[3].y * m[0].z - m[3].x * m[2].y * m[0].z
			+ m[3].x * m[0].y * m[2].z - m[0].x * m[3].y * m[2].z
			- m[2].x * m[0].y * m[3].z + m[0].x * m[2].y * m[3].z;

	r[2].x = m[1].y * m[3].z * m[0].w - m[3].y * m[1].z * m[0].w
			+ m[3].y * m[0].z * m[1].w - m[0].y * m[3].z * m[1].w
			- m[1].y * m[0].z * m[3].w + m[0].y * m[1].z * m[3].w;

	r[2].y = m[3].x * m[1].z * m[0].w - m[1].x * m[3].z * m[0].w
			- m[3].x * m[0].z * m[1].w + m[0].x * m[3].z * m[1].w
			+ m[1].x * m[0].z * m[3].w - m[0].x * m[1].z * m[3].w;

	r[2].z = m[1].x * m[3].y * m[0].w - m[3].x * m[1].y * m[0].w
			+ m[3].x * m[0].y * m[1].w - m[0].x * m[3].y * m[1].w
			- m[1].x * m[0].y * m[3].w + m[0].x * m[1].y * m[3].w;

	r[2].w = m[3].x * m[1].y * m[0].z - m[1].x * m[3].y * m[0].z
			- m[3].x * m[0].y * m[1].z + m[0].x * m[3].y * m[1].z
			+ m[1].x * m[0].y * m[3].z - m[0].x * m[1].y * m[3].z;

	r[3].x = m[2].y * m[1].z * m[0].w - m[1].y * m[2].z * m[0].w
			- m[2].y * m[0].z * m[1].w + m[0].y * m[2].z * m[1].w
			+ m[1].y * m[0].z * m[2].w - m[0].y * m[1].z * m[2].w;

	r[3].y = m[1].x * m[2].z * m[0].w - m[2].x * m[1].z * m[0].w
			+ m[2].x * m[0].z * m[1].w - m[0].x * m[2].z * m[1].w
			- m[1].x * m[0].z * m[2].w + m[0].x * m[1].z * m[2].w;

	r[3].z = m[2].x * m[1].y * m[0].w - m[1].x * m[2].y * m[0].w
			- m[2].x * m[0].y * m[1].w + m[0].x * m[2].y * m[1].w
			+ m[1].x * m[0].y * m[2].w - m[0].x * m[1].y * m[2].w;

	r[3].w = m[1].x * m[2].y * m[0].z - m[2].x * m[1].y * m[0].z
			+ m[2].x * m[0].y * m[1].z - m[0].x * m[2].y * m[1].z
			- m[1].x * m[0].y * m[2].z + m[0].x * m[1].y * m[2].z;

	r *= det;

	return r;
}

bool linePlaneIntersection(Angel::vec3 lineStart, Angel::vec3 lineEnd,
		Angel::vec3 pointOnPlane, Angel::vec3 planeNormal,
		Angel::vec3 &intersection, float &dist) {
	Angel::vec3 lineDir = lineEnd - lineStart;
	lineDir = Angel::normalize(lineDir);

	Angel::vec3 planeDir = pointOnPlane - lineStart;

	float n = Angel::dot(planeNormal, planeDir);
	float d = Angel::dot(planeNormal, lineDir);

	if (fabs(d) < 0.00001) {
		// Line is parallel to plane
		return false;
	}

	dist = n / d;

	Angel::vec3 offset = dist * lineDir;
	intersection = lineStart + offset;

	return true;
}

void projectScreenPointToPlane(Angel::mat4 modelViewMatrix,
		Angel::mat4 inverseProjMatrix, int x, int y, Angel::vec3 planeCenter,
		Angel::vec3 planeNormal, Angel::vec3 &intersection, float &dist) {
	float fx = ((float) x / SCREEN_WIDTH * 2.0f - 1.0);
	float fy = -((float) y / SCREEN_HEIGHT * 2.0f - 1.0);

	Angel::vec4 ndcNear(fx, fy, -1, 1);
	Angel::vec4 ndcFar(fx, fy, 1, 1);

	// Normalized Device Coordinates to Eye Coordinates
	Angel::vec4 pointOnNearPlane = inverseProjMatrix * ndcNear;
	Angel::vec4 pointOnFarPlane = inverseProjMatrix * ndcFar;
	pointOnNearPlane = pointOnNearPlane / pointOnNearPlane.w;
	pointOnFarPlane = pointOnFarPlane / pointOnFarPlane.w;

	// Eye Coordinates to Object Coordinates
	Angel::mat4 inverseModelViewMatrix = inverse4(modelViewMatrix);

	Angel::vec4 nearWorld = inverseModelViewMatrix * pointOnNearPlane;
	Angel::vec4 farWorld = inverseModelViewMatrix * pointOnFarPlane;

	Angel::vec3 lineStart = Angel::vec3(nearWorld.x, nearWorld.y, nearWorld.z);
	Angel::vec3 lineEnd = Angel::vec3(farWorld.x, farWorld.y, farWorld.z);

	linePlaneIntersection(lineStart, lineEnd, planeCenter, planeNormal,
			intersection, dist);
}

//----------------------------------------------------------------------------

// http://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Arcball
/**
 * Get a normalized vector from the center of the virtual ball O to a
 * point P on the virtual ball surface, such that P is aligned on
 * screen's (X,Y) coordinates.  If (X,Y) is too far away from the
 * sphere, return the nearest point on the virtual ball surface.
 */
Angel::vec3 get_arcball_vector(int x, int y) {
	Angel::vec3 P = Angel::vec3(1.0 * x / SCREEN_WIDTH * 2 - 1.0,
			1.0 * y / SCREEN_HEIGHT * 2 - 1.0, 0);
	P.y = -P.y;
	float OP_squared = P.x * P.x + P.y * P.y;
	if (OP_squared <= 1 * 1)
		// Note that P.z is flipped due to left-handedness
		P.z = sqrt(1 * 1 - OP_squared);  // Pythagorean theorem
	else
		P = Angel::normalize(P);  // nearest point
	return P;
}

Angel::mat4 angleAxis(float angle, Angel::vec3 axis) {
	// http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
	float c = cos(angle);
	float s = sin(angle);
	float C = 1 - c;

	return Angel::mat4(
			c + pow(axis.x, 2) * (1-c), axis.x*axis.y*(1-c)-axis.z*s, axis.x*axis.z*(1-c)+axis.y*s, 0.0,
			axis.y*axis.x*(1-c)+axis.z*s, c + pow(axis.y, 2)*(1-c), axis.y*axis.z*(1-c)-axis.x*s, 0.0,
			axis.z*axis.x*(1-c)-axis.y*s, axis.z*axis.y*(1-c)+axis.x*s, c + pow(axis.z, 2)*(1-c), 0.0,
			0.0, 0.0, 0.0, 1.0
		);
}

//----------------------------------------------------------------------------

void copyMat4to16(Angel::mat4 &m, GLfloat a[16]) {
	// Not the best way to copy a 4x4 matrix to a 16-element array
	// but the most straightforward
	a[0] = m[0][0];
	a[1] = m[1][0];
	a[2] = m[2][0];
	a[3] = m[3][0];
	a[4] = m[0][1];
	a[5] = m[1][1];
	a[6] = m[2][1];
	a[7] = m[3][1];
	a[8] = m[0][2];
	a[9] = m[1][2];
	a[10] = m[2][2];
	a[11] = m[3][2];
	a[12] = m[0][3];
	a[13] = m[1][3];
	a[14] = m[2][3];
	a[15] = m[3][3];
}

//----------------------------------------------------------------------------

// quad generates two triangles for each face and assigns colors
//    to the vertices
int Index = 0;
void quad(int a, int b, int c, int d, Angel::mat4 translate, bool black=false) {
	int col;

	if (black)
		col = 0;
	else
		col = a;

	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[a];
	Index++;
	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[b];
	Index++;
	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[c];
	Index++;
	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[a];
	Index++;
	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[c];
	Index++;
	oColors[Index] = vertex_colors[col];
	oPoints[Index] = translate * vertices[d];
	Index++;
}

//----------------------------------------------------------------------------

// generate 12 triangles: 36 vertices and 36 colors
void rubixcube() {
	float locs[3] = {-CUBE_OFFSET, 0, CUBE_OFFSET};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				quad(1, 0, 3, 2, Angel::Translate(locs[i], locs[j], locs[k]), k < 2); // white
				quad(2, 3, 7, 6, Angel::Translate(locs[i], locs[j], locs[k]), i < 2); // red
				quad(3, 0, 4, 7, Angel::Translate(locs[i], locs[j], locs[k]), j > 0); // blue
				quad(6, 5, 1, 2, Angel::Translate(locs[i], locs[j], locs[k]), j < 2); // green
				quad(4, 5, 6, 7, Angel::Translate(locs[i], locs[j], locs[k]), k > 0); // yellow
				quad(5, 4, 0, 1, Angel::Translate(locs[i], locs[j], locs[k]), i > 0); // orange
			}
		}
	}

	std::copy(&oPoints[0], &oPoints[NumVertices], points);
	std::copy(&oColors[0], &oColors[NumVertices], colors);
}

void rotateCube(int index1, int index2, int axis) {
	// rotates the cube specified by index1 and stores it in the location specified by index2
	switch(axis) {
	case Xaxis:
		for (int v = 0; v < 6; v++) {
			colors[index1 + 6*0 + v] = oColors[index2 + 6*1 + v];
			colors[index1 + 6*1 + v] = oColors[index2 + 6*4 + v];
			colors[index1 + 6*4 + v] = oColors[index2 + 6*5 + v];
			colors[index1 + 6*5 + v] = oColors[index2 + 6*0 + v];
			colors[index1 + 6*2 + v] = oColors[index2 + 6*2 + v];
			colors[index1 + 6*3 + v] = oColors[index2 + 6*3 + v];
		}
		break;
	case Yaxis:
		for (int v = 0; v < 6; v++) {
			colors[index1 + 6*0 + v] = oColors[index2 + 6*3 + v];
			colors[index1 + 6*3 + v] = oColors[index2 + 6*4 + v];
			colors[index1 + 6*4 + v] = oColors[index2 + 6*2 + v];
			colors[index1 + 6*2 + v] = oColors[index2 + 6*0 + v];
			colors[index1 + 6*1 + v] = oColors[index2 + 6*1 + v];
			colors[index1 + 6*5 + v] = oColors[index2 + 6*5 + v];
		}
		break;
	case Zaxis:
		for (int v = 0; v < 6; v++) {
			colors[index1 + 6*1 + v] = oColors[index2 + 6*3 + v];
			colors[index1 + 6*3 + v] = oColors[index2 + 6*5 + v];
			colors[index1 + 6*5 + v] = oColors[index2 + 6*2 + v];
			colors[index1 + 6*2 + v] = oColors[index2 + 6*1 + v];
			colors[index1 + 6*4 + v] = oColors[index2 + 6*4 + v];
			colors[index1 + 6*0 + v] = oColors[index2 + 6*0 + v];
		}
		break;
	}
}

bool isColorEqual(Angel::vec4 c1, Angel::vec4 c2) {
	// checks if two colors are equal
	return c1.x == c2.x && c1.y == c2.y && c1.z == c2.z && c1.w == c2.w;
}

bool checkSolved() {
	// iterates through the cube's faces and checkes to make sure that
	// the cublet faces are all of the same color

	for (int axis = 0; axis < 3; axis++) {
		// j
		for (int j = 0; j < 3; j += 2) {
			// only need to check outer faces
			color4 c;
			int face;
			bool assign = false;
			if (j == 0)
				face = 2;
			else
				face = 3;

			for (int i = 0; i < 3; i++) {
				for (int k = 0; k < 3; k++) {
					// note all vertices have the same color
					if (!assign) {
						assign = true;
						c = colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face];
					} else if (!isColorEqual(colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face], c)) {
						return false;
					}
				}
			}
		}
		// i
		for (int i = 0; i < 3; i += 2) {
			// only need to check outer faces

			color4 c;
			int face;
			bool assign = false;
			if (i == 0)
				face = 5;
			else
				face = 1;

			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					// note all vertices have the same color
					if (!assign) {
						assign = true;
						c = colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face];
					} else if (!isColorEqual(colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face], c)) {
						return false;
					}
				}
			}
		}
		// k
		for (int k = 0; k < 3; k += 2) {
			// only need to check outer faces

			color4 c;
			int face;
			bool assign = false;
			if (k == 0)
				face = 4;
			else
				face = 0;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					// note all vertices have the same color
					if (!assign) {
						assign = true;
						c = colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face];
					} else if (!isColorEqual(colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*face], c)) {
						return false;
					}
				}
			}
		}
	}

	return true;
}

//----------------------------------------------------------------------------

// OpenGL initialization
void init() {
	srand(time(0));

	rubixcube();

	// Create a vertex array object
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// Create and initialize a buffer object
	GLuint buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(points) + sizeof(colors), NULL,
			GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(points), points);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(points), sizeof(colors), colors);

	// Load shaders and use the resulting shader program
	GLuint program = InitShader("vshader.glsl", "fshader.glsl");
	glUseProgram(program);

	// set up vertex arrays
	GLuint vPosition = glGetAttribLocation(program, "vPosition");
	glEnableVertexAttribArray(vPosition);
	glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0,
			BUFFER_OFFSET(0));

	GLuint vColor = glGetAttribLocation(program, "vColor");
	glEnableVertexAttribArray(vColor);
	glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0,
			BUFFER_OFFSET(sizeof(points)));

	arcball = glGetUniformLocation(program, "arcball");
	zoom = glGetUniformLocation(program, "zoom");

	copyMat4to16(mZoom, Zoom);
	copyMat4to16(mArcball, Arcball);

	glEnable (GL_DEPTH_TEST);
	glClearColor(.5, .5, .5, 1.0);
}

//----------------------------------------------------------------------------

void display(void) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUniformMatrix4fv(arcball, 1, true, Arcball);
	glUniformMatrix4fv(zoom, 1, true, Zoom);
	glDrawArrays(GL_TRIANGLES, 0, NumVertices);

	glutSwapBuffers();
}

//----------------------------------------------------------------------------

void resetColors() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				// rotates the layer by 90 deg
				for (int c = 0; c < 6; c++) {
					for (int v = 0; v < 6; v++) {
						colors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*c + v] =
								oColors[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*c + v];
					}
				}
			}
		}
	}
}

//----------------------------------------------------------------------------

void rotateLayer(int axis, int layer, int quarters) {
	// rotates the layer around the axis a specified number of 90 degree turns
	for (int r = 0; r < quarters; r++) {
		switch (axis) {
		case Xaxis:
			int j = layer;
			for (int i = 0; i < 3; i++) {
				for (int k = 0; k < 3; k++) {
					// rotates the layer by 90 deg
					rotateCube(6*6*3*3*i + 6*6*3*j + 6*6*k, 6*6*3*3*k + 6*6*3*j + 6*6*abs(i-2), Xaxis);
				}
			}
			break;
		case Yaxis:
			int i = layer;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					// rotates the layer by 90 deg
					rotateCube(6*6*3*3*i + 6*6*3*j + 6*6*k, 6*6*3*3*i + 6*6*3*k + 6*6*abs(j-2), Yaxis);
				}
			}
			break;
		case Zaxis:
			int k = layer;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					// rotates the layer by 90 deg
					rotateCube(6*6*3*3*i + 6*6*3*j + 6*6*k, 6*6*3*3*abs(j-2) + 6*6*3*i + 6*6*k, Zaxis);
				}
			}
			break;
		}

		std::copy(&colors[0], &colors[NumVertices], oColors);
	}

	if (checkSolved()) {
		glClearColor(.5, .5, .5, 1.0);
	} else {
		glClearColor(.9, .9, .9, 1.0);
	}

	glBufferSubData(GL_ARRAY_BUFFER, sizeof(points), sizeof(colors), colors);
}

//----------------------------------------------------------------------------

void twistLayer(int i, int j, int k) {
	// Twists the specified layer a certain number of degrees determined by mTwist
	for (int c = 0; c < 6; c++) {
		for (int v = 0; v < 6; v++) {
			points[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*c + v] =
					mTwist * oPoints[6*6*3*3*i + 6*6*3*j + 6*6*k + 6*c + v];
		}
	}
}

void resetTwist() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				twistLayer(i, j, k);
			}
		}
	}
}

void twistCube(twister t) {
	if (t.axis.x >= 1 || t.axis.x <= -1) {
		int i = t.layer;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				twistLayer(i, j, k);
			}
		}
	} else if (t.axis.y >= 1 || t.axis.y <= -1) {
		// axis is flipped because of left->right handed coordinate system
		int j = t.layer;
		for (int i = 0; i < 3; i++) {
			for (int k = 0; k < 3; k++) {
				twistLayer(i, j, k);
			}
		}
	} else if (t.axis.z >= 1 || t.axis.z <= -1) {
		int k = t.layer;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				twistLayer(i, j, k);
			}
		}
	} else {
		resetTwist();
	}
}

int snap(GLfloat angle) {
	while (angle < 0)
		angle += 2*M_PI;
	while (angle >= 2*M_PI)
		angle -= 2*M_PI;

	if (angle <= M_PI/4 && angle >= -M_PI/4)
		return 0;
	else if (angle > M_PI/4 && angle < 3*M_PI/4)
		return 1;
	else if (angle >= 3*M_PI/4 && angle <= 5*M_PI/4)
		return 2;
	else if (angle > 5*M_PI/4 && angle < 7*M_PI/4)
		return 3;
	else if (angle >= 7*M_PI/4 && angle < 9*M_PI/4)
		return 0;
	else
		return 0;
}

//----------------------------------------------------------------------------

void save(char* file) {
	// saves the current state to the given file

	std::ofstream saveFile;
	saveFile.open(file);

	for (int i = 0; i < NumVertices; i++) {
		saveFile << oColors[i].x << "," << oColors[i].y << "," << oColors[i].z << "," << oColors[i].w << std::endl;
	}

	saveFile.close();

	std::cout << "Saved to: " << file << std::endl;
}

float parseFloat(std::string s) {
	std::istringstream ss(s);
	float f;
	ss >> f;
	return f;
}

void load(char* file) {
	// loads the state from the given file
	std::ifstream loadFile(file);
	if (loadFile.is_open()) {
		std::string line;
		int i = 0;
		while (loadFile.good()) {
			getline(loadFile, line);
			std::istringstream liness(line);

			std::string xs, ys, zs, ws;
			getline(liness, xs, ',');
			getline(liness, ys, ',');
			getline(liness, zs, ',');
			getline(liness, ws, ',');

			oColors[i] = point4(parseFloat(xs), parseFloat(ys), parseFloat(zs), parseFloat(ws));
			i++;
		}

		loadAFile = false;
		loadFile.close();

		std::cout << "Loaded: " << file << std::endl;
	}
	else
		std::cout << "Unable to open file" << std::endl;

	// display
	std::copy(&oColors[0], &oColors[NumVertices], colors);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(points), sizeof(colors), colors);

	if (checkSolved()) {
		glClearColor(.5, .5, .5, 1.0);
	} else {
		glClearColor(.9, .9, .9, 1.0);
	}

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 033: // Escape Key
	case 'q':
	case 'Q':
		exit (EXIT_SUCCESS);
		break;
	case ' ':
		mArcball = RotateX(45) * RotateY(45) * RotateZ(45);
		copyMat4to16(mArcball, Arcball);
		glutPostRedisplay();
		break;
	case 's':
		save(saveName);
		break;
	case 13:
		rotateLayer(rand() % 3, rand() % 3, rand() % 3 + 1);
		glutPostRedisplay();
		break;
	}
}

//----------------------------------------------------------------------------

void mouse(int button, int state, int x, int y) {
	if (state == GLUT_DOWN) {
		switch (button) {
		case GLUT_LEFT_BUTTON:
			if (!mouse2) {
				mouse1 = true;
				oldX = x;
				oldY = y;
			}
			break;
		case GLUT_RIGHT_BUTTON:
			if (!mouse1 && !mouse2) {
				mouse2 = true;
				twist.angle = 0;
				twist.dir = -1;

				plane select;
				Angel::vec3 intersect, closest;
				float d;
				float shortest = FLT_MIN;

				for (int i = 0; i < 3; i++) {

					Angel::mat4 rotate;

					// rotate to view different perspectives
					switch(i) {
					case Xaxis: rotate = Angel::RotateX(90) * Angel::Scale(-1, -1, 1); // scale switches to left-handed coordinates
					break;
					case Yaxis: rotate = Angel::RotateY(90) * Angel::Scale(-1, -1, 1); // scale switches to left-handed coordinates
					break;
					case Zaxis: rotate = Angel::mat4();
					break;
					}

					// iterate through front and back
					for (int j = -1; j <= 1; j += 2) {
						Angel::vec3 axis = Angel::vec3(0, 0, 1) * (float) j;

						plane current;
						current.normal = axis;
						current.point = axis * CUBE_SIZE / 2;
						current.rotate = rotate;
						current.rotation = i;

						projectScreenPointToPlane(current.rotate * mZoom * mArcball, Angel::Ortho(-1, 1, -1, 1, -1, 1), x, y,
								current.point, current.normal, intersect, d);

						if (intersect.x >= -CUBE_SIZE/2 && intersect.x <= CUBE_SIZE/2 &&
								intersect.y >= -CUBE_SIZE/2 && intersect.y <= CUBE_SIZE/2) {
							if (shortest == FLT_MIN || d > shortest) {
								closest = intersect;
								shortest = d;
								select = current;
							}
						}
					}
				}

				if (shortest != FLT_MIN) {
					selected = true;
					selectedPlane = select;
					selectStart = closest;
				}
			}
			break;
		}
	}

	if (state == GLUT_UP) {
		switch (button) {
		case GLUT_LEFT_BUTTON:
			mouse1 = false;
			break;
		case GLUT_RIGHT_BUTTON:
			mouse2 = false;

			if (selected && twist.layer != -1) {
				// adjust angle for axis
				int s = sign(round(twist.axis.x)) - sign(round(twist.axis.y)) - sign(round(twist.axis.z));
				// snap the movement to a 90 degree rotation
				int quarters = snap(s * twist.angle);

				if (quarters != 0) {
					int axis = abs((int) round(twist.axis.x)) + 2*abs((int) round(twist.axis.z));
					rotateLayer(axis, twist.layer, quarters);
				}
			}

			mTwist = Angel::mat4();
			resetTwist();
			selected = false;

			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(points), points);
			glutPostRedisplay();

			break;
		}
	}
}

void motion(int x, int y) {
	if (mouse1) {
		// http://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Arcball
		if (x != oldX || y != oldY) {
			Angel::vec3 va = get_arcball_vector(oldX, oldY);
			Angel::vec3 vb = get_arcball_vector(x, y);

			Angel::vec3 axis = Angel::normalize(Angel::cross(va, vb));
			float angle = std::acos(std::min(1.0f, Angel::dot(va, vb)));

			// save old orientation
			mat4 oldArcball = mArcball;

			// convert axis angle representation to angles
			mArcball = angleAxis(angle, axis);

			// TODO: a more complete check
			if (isNumber(mArcball[0].x)) { // check to ensure matrix is valid
				// rotate based on previous rotation
				mArcball = oldArcball * mArcball;

				copyMat4to16(mArcball, Arcball);

				oldX = x;
				oldY = y;
			} else {
				// pass through
				mArcball = oldArcball;
			}
		}

		glutPostRedisplay();
	} else if (mouse2 && selected) {
		// We need to determine which of the two orthogonal planes
		// (to the selected plane) the user wants to rotate on.

		// This is done by projecting the movement vector onto the plane
		// and selecting the orthogonal axis the projection is closer to.

		Angel::vec3 selectEnd, direction;
		float dist;
		int layer = -1;

		projectScreenPointToPlane(selectedPlane.rotate * mZoom * mArcball, Angel::Ortho(-1, 1, -1, 1, -1, 1), x, y,
				selectedPlane.point, selectedPlane.normal, selectEnd, dist);

		Angel::vec3 disp = selectEnd - selectStart;
		disp.z = 0; // it should be zero

		if (Angel::length(disp) >= TWIST_MIN_DIST) {

			Angel::vec3 axis1 = orthoShift(selectedPlane.normal);
			Angel::vec3 axis2 = orthoShift(axis1);
			Angel::vec3 axis, layerAxis;

			float proj1 = Angel::dot(Angel::normalize(disp), Angel::normalize(axis1));
			float proj2 = Angel::dot(Angel::normalize(disp), Angel::normalize(axis2));

			float angle1 = acos(proj1) * 180/M_PI;
			float angle2 = acos(proj2) * 180/M_PI;

			// screen distance determines angle
			float dist = Angel::length(Angel::vec2(x, y) - Angel::vec2(oldX, oldY));

			// pick the correct axis
			if ((twist.dir == -1 || twist.dir == 0) &&
					(angle1 < TWIST_MAX_ANGLE || angle1 > 180 - TWIST_MAX_ANGLE)) {
				axis = axis1;
				layerAxis = axis2;

				// little bit of a hack to account for certain faces
				float proj = -(selectedPlane.rotation == 0 || selectedPlane.rotation == 1 ? -1 : 1) *
						selectStart.y / CUBE_SIZE + .5;
				if (proj < 1.0f/3)
					layer = 2;
				else if (proj >= 1.0f/3 && proj <= 2.0f/3)
					layer = 1;
				else if (proj > 2.0f/3)
					layer = 0;
				if (layer != -1)
					twist.dir = 0;
			} else if ((twist.dir == -1 || twist.dir == 1) &&
					(angle2 < TWIST_MAX_ANGLE || angle2 > 180 - TWIST_MAX_ANGLE)) {
				axis = axis2;
				layerAxis = axis1;

				float proj = -(selectedPlane.rotation == 0 ? -1 : 1) *
						selectStart.x / CUBE_SIZE + .5;
				if (proj < 1.0f/3)
					layer = 2;
				else if (proj >= 1.0f/3 && proj <= 2.0f/3)
					layer = 1;
				else if (proj > 2.0f/3)
					layer = 0;
				if (layer != -1)
					twist.dir = 1;
			}

			if (layer != -1) {
				direction = transform3(selectedPlane.rotate, Angel::cross(selectedPlane.normal, axis));

				// save axis-angle
				GLfloat angle = Angel::dot(disp, axis);
				// convert axis angle representation to rotation matrix
				Angel::mat4 tTwist = angleAxis(angle, direction);

				// TODO: a more complete check
				if (isNumber(mTwist[0].x)) { // check to ensure matrix is valid
					mTwist = tTwist;
					twist.angle = angle;
					twist.axis = direction;
					twist.layer = layer;

					twistCube(twist);
				}
			} else {
				mTwist = Angel::mat4();
				twist.layer = -1;
				twist.dir = -1;
				resetTwist();
			}
		} else {
			mTwist = Angel::mat4();
			twist.layer = -1;
			twist.dir = -1;
			resetTwist();
		}

		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(points), points);
		glutPostRedisplay();
	}
}

void mouseWheel(int wheel, int direction, int x, int y) {
	if (direction < 0) {
		scale *= 1.0f/SCALE_FACTOR;
		if (scale < 1.0f/12)
			scale = 1.0f/12;
	} else if (direction > 0) {
		scale *= SCALE_FACTOR;
		if (scale > 1.0f/3)
			scale = 1.0f/3;
	}
	mZoom = Angel::Scale(scale, scale, scale);
	copyMat4to16(mZoom, Zoom);
	glutPostRedisplay();
}

//----------------------------------------------------------------------------

void idle(void) {
	if (turns > 0) {
		for (int i = 0; i < turns; i++) {
			rotateLayer(rand() % 3, rand() % 3, rand() % 3 + 1);
		}

		turns = 0;
	}

	if (loadAFile) {
		load(loadName);
	}

	glutPostRedisplay();
}

//----------------------------------------------------------------------------

int main(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_HEIGHT, SCREEN_WIDTH);
//    glutInitContextVersion( 3, 2 );
//    glutInitContextProfile( GLUT_CORE_PROFILE );
	glutCreateWindow("Rubik's Cube");

	glewInit();

	init();

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMouseWheelFunc(mouseWheel);
	glutIdleFunc(idle);

	if (argc > 1) {
		saveName = argv[1];
	} else {
		saveName = "save.txt";
	}

	if (argc > 2) {
		std::istringstream ss(argv[2]);
		int num;

		if (!(ss >> num)) {
			loadName = argv[2];
			loadAFile = true;
		} else {
			turns = num;
		}
	}

	glutMainLoop();

	return 0;
}
