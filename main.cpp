// OpenCSG - library for image-based CSG rendering for OpenGL
// Copyright (C) 2002-2014, Florian Kirsch,
// Hasso-Plattner-Institute at the University of Potsdam, Germany
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License,
// Version 2, as published by the Free Software Foundation.
// As a special exception, you have permission to link this library
// with the CGAL library and distribute executables.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

//
// main.cpp
//
// simple example program for OpenCSG using Glut
//

#include <GL/glew.h>
#include <opencsg.h>
#include "displaylistPrimitive.h"
#include "vboPrimitive.h"
#include <iostream>
#include <sstream>
#include <string>

// include glut.h after stdlib.h to avoid conflict in declaration
// of exit() with Visual Studio 2010
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

enum {
    CSG_BASIC, CSG_WIDGET, CSG_GRID2D, CSG_GRID3D, CSG_CUBERACK, CSG_CONCAVE,

    ALGO_AUTOMATIC, GF_STANDARD, GF_DC, GF_OQ, SCS_STANDARD, SCS_DC, SCS_OQ,

    OFFSCREEN_AUTOMATIC, OFFSCREEN_FBO, OFFSCREEN_PBUFFER,

    USE_DISPLAY_LISTS, USE_VAO_VBOS
};

void checkGLErrors(std::string location)
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR) {
		std::cerr << "OpenGL Error: " << location << std::endl << gluErrorString(err) << std::endl;
	}
}

std::vector<OpenCSG::Primitive*> primitives;

bool		   useVaos = false;
bool		   use_display_lists = true;
bool               spin = true;
float              rot = 0.0f;
std::ostringstream fpsStream;
GLuint faceshader_prog;
std::map<std::string, GLuint> shader_attributes;

// storage for Matrices
std::vector<float> proj_matrix(16);
std::vector<float> view_matrix(16);

void init();

// ----------------------------------------------------
// VECTOR STUFF
//

// res = a cross b;
void crossProduct(float *a, float *b, float *res)
{

	res[0] = a[1] * b[2]  -  b[1] * a[2];
	res[1] = a[2] * b[0]  -  b[2] * a[0];
	res[2] = a[0] * b[1]  -  b[0] * a[1];
}

// Normalize a vec3
void normalize(float *a)
{

	float mag = sqrt(a[0] * a[0]  +  a[1] * a[1]  +  a[2] * a[2]);

	a[0] /= mag;
	a[1] /= mag;
	a[2] /= mag;
}

// ----------------------------------------------------
// MATRIX STUFF
//

// sets the square matrix mat to the identity matrix,
// size refers to the number of rows (or columns)
void setIdentityMatrix(std::vector<float> &mat, int size)
{

	// fill matrix with 0s
	for (int i = 0; i < size * size; ++i)
		mat[i] = 0.0f;

	// fill diagonal with 1s
	for (int i = 0; i < size; ++i)
		mat[i + i * size] = 1.0f;
}

//
// a = a * b;
//
void multMatrix(std::vector<float> &a, std::vector<float> &b)
{

	std::vector<float> res(16, 0);

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			res[j*4 + i] = 0.0f;
			for (int k = 0; k < 4; ++k) {
				res[j*4 + i] += a[k*4 + i] * b[j*4 + k];
			}
		}
	}
	a = res;
}

// Defines a transformation matrix mat with a translation
void setTranslationMatrix(std::vector<float> &mat, float x, float y, float z)
{

	setIdentityMatrix(mat,4);
	mat[12] = x;
	mat[13] = y;
	mat[14] = z;
}

// ----------------------------------------------------
// Projection Matrix
//

void buildProjectionMatrix(float fov, float ratio, float nearP, float farP)
{

	float f = 1.0f / tan (fov * (M_PI / 360.0));

	setIdentityMatrix(proj_matrix,4);

	proj_matrix[0] = f / ratio;
	proj_matrix[1 * 4 + 1] = f;
	proj_matrix[2 * 4 + 2] = (farP + nearP) / (nearP - farP);
	proj_matrix[3 * 4 + 2] = (2.0f * farP * nearP) / (nearP - farP);
	proj_matrix[2 * 4 + 3] = -1.0f;
	proj_matrix[3 * 4 + 3] = 0.0f;
}

// ----------------------------------------------------
// View Matrix
//
// note: it assumes the camera is not tilted,
// i.e. a vertical up vector (remmeber gluLookAt?)
//

void setCamera(float posX, float posY, float posZ, float lookAtX, float lookAtY, float lookAtZ)
{

	float dir[3], right[3], up[3];

	up[0] = 0.0f;	up[1] = 1.0f;	up[2] = 0.0f;

	dir[0] =  (lookAtX - posX);
	dir[1] =  (lookAtY - posY);
	dir[2] =  (lookAtZ - posZ);
	normalize(dir);

	crossProduct(dir,up,right);
	normalize(right);

	crossProduct(right,dir,up);
	normalize(up);

	std::vector<float> aux(16, 0);

	view_matrix[0]  = right[0];
	view_matrix[4]  = right[1];
	view_matrix[8]  = right[2];
	view_matrix[12] = 0.0f;

	view_matrix[1]  = up[0];
	view_matrix[5]  = up[1];
	view_matrix[9]  = up[2];
	view_matrix[13] = 0.0f;

	view_matrix[2]  = -dir[0];
	view_matrix[6]  = -dir[1];
	view_matrix[10] = -dir[2];
	view_matrix[14] =  0.0f;

	view_matrix[3]  = 0.0f;
	view_matrix[7]  = 0.0f;
	view_matrix[11] = 0.0f;
	view_matrix[15] = 1.0f;

	setTranslationMatrix(aux, -posX, -posY, -posZ);

	multMatrix(view_matrix, aux);
}

void changeSize(int w, int h)
{

	float ratio;
	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	ratio = (1.0f * w) / h;
	buildProjectionMatrix(53.13f, ratio, 1.0f, 30.0f);
}


/////////////////////

void clearPrimitives()
{
	for (std::vector<OpenCSG::Primitive*>::const_iterator i = primitives.begin(); i != primitives.end(); ++i) {
		OpenCSG::Primitive* p =
			static_cast<OpenCSG::Primitive*>(*i);
		delete p;
	}

	primitives.clear();
}

void solidCylinder(GLdouble radius, GLdouble height, GLint slices, GLint stacks)
{

	GLUquadricObj* qobj = gluNewQuadric();

	gluCylinder(qobj, radius, radius, height, slices, stacks);
	glScalef(-1.0f, 1.0f, -1.0f);
	gluDisk(qobj, 0.0, radius, slices, stacks);
	glScalef(-1.0f, 1.0f, -1.0f);
	glTranslatef(0.0f, 0.0f, static_cast<GLfloat>(height));
	gluDisk(qobj, 0.0, radius, slices, stacks);

	gluDeleteQuadric(qobj);
}

void setBasicShape()
{
	clearPrimitives();

	if (use_display_lists) {
		GLuint id1 = glGenLists(1);
		glNewList(id1, GL_COMPILE);
		glPushMatrix();
		glTranslatef(-0.25f, 0.0f, 0.0f);
		glutSolidSphere(1.0, 20, 20);
		glPopMatrix();
		glEndList();

		GLuint id2 = glGenLists(1);
		glNewList(id2, GL_COMPILE);
		glPushMatrix();
		glTranslatef(0.25f, 0.0f, 0.0f);
		glutSolidSphere(1.0, 20, 20);
		glPopMatrix();
		glEndList();

		GLuint id3 = glGenLists(1);
		glNewList(id3, GL_COMPILE);
		glPushMatrix();
		glTranslatef(0.0f, 0.0f, 0.5f);
		glScalef(0.5f, 0.5f, 2.0f);
		glutSolidSphere(1.0, 20, 20);
		glPopMatrix();
		glEndList();

		primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id2, OpenCSG::Intersection, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id3, OpenCSG::Subtraction, 1));
	} else {
		primitives.push_back(new OpenCSG::SpherePrimitive(1.0, false, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate(-0.25f, 0.0f, 0.0f);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::SpherePrimitive(1.0, false, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate(0.25f, 0.0f, 0.0f);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::SpherePrimitive(1.0, false, OpenCSG::Subtraction, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate(0.0f, 0.0f, 0.5f);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->scale(0.5f, 0.5f, 2.0f);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
	}
}

void setWidget()
{
	clearPrimitives();

	if (use_display_lists) {
		GLuint id1 = glGenLists(1);
		glNewList(id1, GL_COMPILE);
		glutSolidSphere(1.2, 20, 20);
		glEndList();

		GLuint id2 = glGenLists(1);
		glNewList(id2, GL_COMPILE);
		glutSolidCube(1.8);
		glEndList();

		GLuint id3 = glGenLists(1);
		glNewList(id3, GL_COMPILE);
		glPushMatrix();
		glTranslatef(0.0f, 0.0f, -1.25f);
		solidCylinder(0.6, 2.5, 20, 20);
		glPopMatrix();
		glEndList();

		GLuint id4 = glGenLists(1);
		glNewList(id4, GL_COMPILE);
		glPushMatrix();
		glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
		glTranslatef(0.0f, 0.0f, -1.25f);
		solidCylinder(0.6, 2.5, 20, 20);
		glPopMatrix();
		glEndList();

		GLuint id5 = glGenLists(1);
		glNewList(id5, GL_COMPILE);
		glPushMatrix();
		glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
		glTranslatef(0.0f, 0.0f, -1.25f);
		solidCylinder(0.6, 2.5, 20, 20);
		glPopMatrix();
		glEndList();

		primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id2, OpenCSG::Intersection, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id3, OpenCSG::Subtraction, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id4, OpenCSG::Subtraction, 1));
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id5, OpenCSG::Subtraction, 1));
	} else {
		primitives.push_back(new OpenCSG::SpherePrimitive(1.2, true, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::CubePrimitive(1.8, 1.8, 1.8, true, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::CylinderPrimitive(0.6, 0.6, 2.5, true, OpenCSG::Subtraction, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::CylinderPrimitive(0.6, 0.6, 2.5, true, OpenCSG::Subtraction, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->rotate(0,90,0);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();

		primitives.push_back(new OpenCSG::CylinderPrimitive(0.6, 0.6, 2.5, true, OpenCSG::Subtraction, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->rotate(90,0,0);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
	}
}

void setGrid2D()
{

	clearPrimitives();

	if (use_display_lists) {
		GLuint id1 = glGenLists(1);
		glNewList(id1, GL_COMPILE);
		glPushMatrix();
		glScalef(1.0f, 0.2f, 1.0f);
		glTranslatef(0.0f, -1.25f, 0.0f);
		glutSolidCube(2.5);
		glPopMatrix();
		glEndList();

		primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
	} else {
		primitives.push_back(new OpenCSG::CubePrimitive(2.5, 2.5, 2.5, true, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->scale(1.0f, 0.2f, 1.0f);
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
	}

	for (int x=-2; x<=2; ++x) {
		for (int z=-2; z<=2; ++z) {
			if (use_display_lists) {
				GLuint id = glGenLists(1);
				glNewList(id, GL_COMPILE);
				glPushMatrix();
				glTranslatef(x*0.5f, 0.0f, z*0.5f);
				glutSolidSphere(0.22, 15, 15);
				glPopMatrix();
				glEndList();

				primitives.push_back(new OpenCSG::DisplayListPrimitive(id, OpenCSG::Subtraction, 1));
			} else {
				primitives.push_back(new OpenCSG::SpherePrimitive(0.22, true, OpenCSG::Subtraction, 1, useVaos));
				dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate(x*0.5f, 0.22f/2.0f+(2.5f*0.2f)/2.0f, z*0.5f);
				dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
			}
		}
	}
}

void setGrid3D()
{
	clearPrimitives();

	if (use_display_lists) {
		GLuint id1 = glGenLists(1);
		glNewList(id1, GL_COMPILE);
		glutSolidCube(2.0);
		glEndList();

		primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
	} else {
		primitives.push_back(new OpenCSG::CubePrimitive(2.0, 2.0, 2.0, true, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
	}

	for (int x=-1; x<=1; ++x) {
		for (int y=-1; y<=1; ++y) {
			for (int z=-1; z<=1; ++z) {
				if (use_display_lists) {
					GLuint id = glGenLists(1);
					glNewList(id, GL_COMPILE);
					glPushMatrix();
					glTranslatef(static_cast<GLfloat>(x), static_cast<GLfloat>(y), static_cast<GLfloat>(z));
					glutSolidSphere(0.58, 20, 20);
					glPopMatrix();
					glEndList();

					primitives.push_back(new OpenCSG::DisplayListPrimitive(id, OpenCSG::Subtraction, 1));
				} else {
					primitives.push_back(new OpenCSG::SpherePrimitive(0.58, false, OpenCSG::Subtraction, 1, useVaos));
					dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate(x,y,z);
					dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
				}
			}
		}
	}
}

void setCubeRack()
{

	clearPrimitives();

	if (use_display_lists) {
		GLuint id1 = glGenLists(1);
		glNewList(id1, GL_COMPILE);
		glutSolidCube(2.0);
		glEndList();

		primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
	} else {
		primitives.push_back(new OpenCSG::CubePrimitive(2.0, 2.0, 2.0, true, OpenCSG::Intersection, 1, useVaos));
		dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
	}

	// mx*x / my*y / mz*z loop all numbers in [-3, 3] in the following order:
	// 3, -3, 2, -2, 1, -1, 0. Compared to the trivial ordering, this makes
	// the CSG rendering less depending on the camera orientation.
	for (int x=3; x>=0; --x) {
		for (int y=3; y>=0; --y) {
			for (int z=3; z>=0; --z) {
				for (int mx=-1; mx<=1 && mx<=x; mx+=2) {
					for (int my=-1; my<=1 && my<=y; my+=2) {
						for (int mz=-1; mz<=1 && mz<=z; mz+=2) {
							if (use_display_lists) {
								GLuint id = glGenLists(1);
								glNewList(id, GL_COMPILE);
								glPushMatrix();
								glTranslatef(float(x*mx)/6.0f, float(y*my)/6.0f, float(z*mz)/6.0f);
								glutSolidSphere(0.58, 20, 20);
								glPopMatrix();
								glEndList();

								primitives.push_back(new OpenCSG::DisplayListPrimitive(id, OpenCSG::Subtraction, 1));
							} else {
								primitives.push_back(new OpenCSG::SpherePrimitive(0.58, false, OpenCSG::Subtraction, 1, useVaos));
								dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->translate((x*mx)/6.0f,(y*my)/6.0f,(z*mz)/6.0f);
								dynamic_cast<OpenCSG::VboPrimitive *>(primitives.back())->initialize();
							}
						}
					}
				}
			}
		}
	}
}

void setConcave()
{

	clearPrimitives();

	GLuint id1 = glGenLists(1);
	glNewList(id1, GL_COMPILE);
	glutSolidTorus(0.6, 1.0, 25, 25);
	glEndList();
	primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 2));

	for (unsigned int i=0; i<4; ++i) {
		GLuint id = glGenLists(1);
		glNewList(id, GL_COMPILE);
		glPushMatrix();
		glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
		glRotatef(i*90.0f + 45.0f, 1.0f, 0.0f, 0.0f);
		glTranslatef(0.0f, 1.0f, 0.0f);
		glutSolidTorus(0.3, 0.6, 15, 15);
		glPopMatrix();
		glEndList();
		primitives.push_back(new OpenCSG::DisplayListPrimitive(id, OpenCSG::Subtraction, 2));
	}

	GLuint id3 = glGenLists(1);
	glNewList(id3, GL_COMPILE);
	glPushMatrix();
	glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
	glTranslatef(0.0f, 0.0f, -1.65f);
	solidCylinder(0.3, 3.3, 20, 20);
	glPopMatrix();
	glEndList();
	primitives.push_back(new OpenCSG::DisplayListPrimitive(id3, OpenCSG::Subtraction, 1));

	GLuint id4 = glGenLists(1);
	glNewList(id4, GL_COMPILE);
	glPushMatrix();
	glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
	glTranslatef(0.0f, 0.0f, -1.65f);
	solidCylinder(0.3, 3.3, 20, 20);
	glPopMatrix();
	glEndList();
	primitives.push_back(new OpenCSG::DisplayListPrimitive(id4, OpenCSG::Subtraction, 1));
}

void renderfps()
{
	glDisable(GL_DEPTH_TEST);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glColor3f(0.0f, 0.0f, 0.0f);
	glRasterPos2f(-1.0f, -1.0f);
	glDisable(GL_LIGHTING);
	std::string s = fpsStream.str();
	for (unsigned int i=0; i<s.size(); ++i) {
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
	}
	glEnable(GL_LIGHTING);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
}

void display()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 2.0, 5.0,  /* eye is at (0,2,5) */
		  0.0, 0.0, 0.0,  /* center is at (0,0,0) */
		  0.0, 1.0, 0.0); /* up is in positive Y direction */
	glRotatef(rot, 0.0f, 1.0f, 0.0f);

	setCamera(0.0, 2.0, 5.0,
		  0.0, 0.0, 0.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	OpenCSG::render(primitives);
	glDepthFunc(GL_EQUAL);
	for (std::vector<OpenCSG::Primitive*>::const_iterator i = primitives.begin(); i != primitives.end(); ++i) {
		(*i)->render();
	}
	glDepthFunc(GL_LESS);

	renderfps();

	glutSwapBuffers();
}

void idle()
{
	static int ancient = 0;
	static int last = 0;
	static int msec = 0;
	last = msec;
	msec = glutGet(GLUT_ELAPSED_TIME);
	if (spin) {
		rot += (msec-last)/10.0f;
		while (rot >= 360.0f)
			rot -= 360.0f;
	}

	static int fps = 0;
	if (last / 1000 != msec / 1000) {

		float correctedFps = static_cast<float>(fps) * 1000.0f / static_cast<float>(msec - ancient);
		fpsStream.str("");
		fpsStream << "fps: " << correctedFps << std::ends;

		ancient = msec;
		fps = 0;
	}

	display();

	++fps;
}

void key(unsigned char k, int, int)
{
	switch (k) {
		case ' ':
			spin = !spin;
			break;
		default:
			break;
	}
	display();
}

void menu(int value)
{
	switch (value) {
		case CSG_BASIC:      setBasicShape();    break;
		case CSG_WIDGET:     setWidget();        break;
		case CSG_GRID2D:     setGrid2D();        break;
		case CSG_GRID3D:     setGrid3D();        break;
		case CSG_CUBERACK:   setCubeRack();      break;
		case CSG_CONCAVE:    setConcave();       break;

		case ALGO_AUTOMATIC: OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Automatic);
			break;
		case GF_STANDARD:    OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::NoDepthComplexitySampling);
			break;
		case GF_DC:          OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::DepthComplexitySampling);
			break;
		case GF_OQ:          OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::OcclusionQuery);
			break;
		case SCS_STANDARD:   OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::NoDepthComplexitySampling);
			break;
		case SCS_DC:         OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::DepthComplexitySampling);
			break;
		case SCS_OQ:         OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
			OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::OcclusionQuery);
			break;

		case OFFSCREEN_AUTOMATIC: OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::AutomaticOffscreenType);
			break;
		case OFFSCREEN_FBO:       OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::FrameBufferObject);
			break;
		case OFFSCREEN_PBUFFER:   OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::PBuffer);
			break;

		case USE_DISPLAY_LISTS:	use_display_lists = true;
			setWidget();
			break;
		case USE_VAO_VBOS: use_display_lists = false;
			setWidget();
			break;
		default: break;
	}
	display();
}

void printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
			std::cerr << infoLog << std::endl;
		free(infoLog);
	}
}

void printProgramInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
			std::cerr << infoLog << std::endl;
		free(infoLog);
	}
}

void init()
{
	GLuint err;
	// gray background
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

	// Enable two OpenGL lights
	GLfloat light_diffuse[]   = { 1.0f,  0.0f,  0.0f,  1.0f};  // Red diffuse light
	GLfloat light_position0[] = {-1.0f, -1.0f, -1.0f,  0.0f};  // Infinite light location
	GLfloat light_position1[] = { 1.0f,  1.0f,  1.0f,  0.0f};  // Infinite light location

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	// Use depth buffering for hidden surface elimination
	glEnable(GL_DEPTH_TEST);

	// Setup the view of the CSG shape
	glMatrixMode(GL_PROJECTION);
	gluPerspective(40.0, 1.0, 1.0, 10.0);

	glMatrixMode(GL_MODELVIEW);

	// starting CSG shape
	setWidget();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);

	glutInitWindowSize(512, 512);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
	glutCreateWindow("OpenCSG example application");

	int err = glewInit();
	if (GLEW_OK != err) {
		// problem: glewInit failed, something is seriously wrong
		std::cerr << "GLEW Error: " << glewGetErrorString(err) << std::endl;
		return 1;
	}

	GLuint major, minor;
	glGetIntegerv(GL_MAJOR_VERSION, (GLint *)&major);
	glGetIntegerv(GL_MINOR_VERSION, (GLint *)&minor);

	if (glewIsSupported("GL_VERSION_3_1") && major >= 3 && minor >= 1) {
		useVaos = true;
	}

	std::cerr << "OpenGL Version: " << major << ", " << minor << std::endl;

	int menuShape     = glutCreateMenu(menu);
	glutAddMenuEntry("Simple",   CSG_BASIC);
	glutAddMenuEntry("Widget",   CSG_WIDGET);
	glutAddMenuEntry("2D-Grid",  CSG_GRID2D);
	glutAddMenuEntry("3D-Grid",  CSG_GRID3D);
	glutAddMenuEntry("Cuberack", CSG_CUBERACK);
	glutAddMenuEntry("Concave (Display List Only)",  CSG_CONCAVE);

	int menuAlgorithm = glutCreateMenu(menu);
	glutAddMenuEntry("Automatic", ALGO_AUTOMATIC);
	glutAddMenuEntry("Goldfeather standard",GF_STANDARD);
	glutAddMenuEntry("Goldfeather depth complexity sampling", GF_DC);
	glutAddMenuEntry("Goldfeather occlusion query", GF_OQ);
	glutAddMenuEntry("SCS standard", SCS_STANDARD);
	glutAddMenuEntry("SCS depth complexity sampling", SCS_DC);
	glutAddMenuEntry("SCS occlusion query", SCS_OQ);

	int menuSettings = glutCreateMenu(menu);
	glutAddMenuEntry("Automatic", OFFSCREEN_AUTOMATIC);
	glutAddMenuEntry("Frame buffer object", OFFSCREEN_FBO);
	glutAddMenuEntry("PBuffer", OFFSCREEN_PBUFFER);

	int menuMode = glutCreateMenu(menu);
	glutAddMenuEntry("Use Display Lists", USE_DISPLAY_LISTS);
	glutAddMenuEntry("Use VAO/VBOs", USE_VAO_VBOS);

	glutCreateMenu(menu);
	glutAddSubMenu("CSG Shapes", menuShape);
	glutAddSubMenu("CSG Algorithms", menuAlgorithm);
	glutAddSubMenu("Settings", menuSettings);
	glutAddSubMenu("Mode", menuMode);

	// connect to right mouse button
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutDisplayFunc(display);
	glutReshapeFunc(changeSize);
	glutKeyboardFunc(key);

	menu(OFFSCREEN_AUTOMATIC);

	glutIdleFunc(idle);
	init();

	try {
		glutMainLoop();
	} catch (std::exception &e) {
		std::cerr << e.what() << std::endl;
	}

	return 0;
}

