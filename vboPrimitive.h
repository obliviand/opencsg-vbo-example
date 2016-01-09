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
// displaylistPrimitive.h
//
// example for a primitive which renders itself using a display list
//

#include <iostream>
#include <map>
#include <cmath>
#include <opencsg.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3d)
using Eigen::Vector3d;
typedef Eigen::Affine3d Transform3d;

const double GRID_FINE   = 0.00000095367431640625;
#define CSGMODE_DIFFERENCE_FLAG 0x10

namespace OpenCSG {
	class VboPrimitive : public Primitive {
		public:
			typedef struct {
				GLfloat x, y, z;
				GLfloat nx, ny, nz;
				GLfloat r, g, b;
			} Vertex3d;

			typedef struct {
				double x, y;
			} Point2d;

			typedef std::vector<Vertex3d> Polygon;
			typedef std::vector<Polygon> Polygons;

			/// An object of this class contains the OpenGL id of a display
			/// list that is compiled by the application. render() just invokes
			/// this display list.
			/// Operation and convexity are just forwarded to the base Primitive class.
			VboPrimitive(Operation operation, unsigned int convexity, bool use_vaos);
			virtual ~VboPrimitive();

			/// Returns the vao id
			GLuint getVaoId() const;

			size_t getTriangleCount() const { return triangle_count_; }

			/// Initialize the vbo
			virtual void initialize() = 0;

			void rotate(float x, float y, float z);
			void translate(float x, float y, float z);
			void scale(float x, float y, float z);

			/// Calls the renderer.
			virtual void render();

		protected:
			void render_surface(const Polygons &polygons);
			// Returns the number of subdivision of a whole circle, given radius and
			// the three special variaſes $fn, $fs and $fa
			inline int get_fragments_from_r(double r, double fn = 150.0, double fs = 2.0, double fa = 12.0)
			{
				// FIXME: It would be better to refuse to create an object. Let's do more strict error handling
				// in future versions of OpenSCAD
				if (r < GRID_FINE || std::isinf(fn) || std::isnan(fn)) {
					return 3;
				}
				if (fn > 0.0) {
					return (int)(fn >= 3 ? fn : 3);
				}
				return (int)ceil(fmax(fmin(360.0 / fa, r*2*M_PI / fs), 5));
			}
			inline void generate_circle(Point2d *circle, double r, int fragments)
			{
				for (int i=0; i<fragments; i++) {
					double phi = (M_PI*2*i) / fragments;
					circle[i].x = r*cos(phi);
					circle[i].y = r*sin(phi);
				}
			}

			inline void checkGLErrors(std::string location)
			{
				GLenum err;

				while ((err = glGetError()) != GL_NO_ERROR) {
					std::cerr << "OpenGL Error: " << location << ": " << gluErrorString(err) << std::endl;
				}
			}


		private:
			GLuint vao_id_;
			GLuint vbo_id_;

			size_t triangle_count_;
			bool use_vaos_;

			Transform3d m_;

			void draw_triangle(Vertex3d *vertices,
					   const Vertex3d &p0, const Vertex3d &p1, const Vertex3d &p2,
					   double e0f, double e1f, double e2f,
					   double nx, double ny, double nz, double nl,
					   double z, bool mirror);
			void gl_draw_triangle(Vertex3d *vertices, const Vertex3d &p0, const Vertex3d &p1, const Vertex3d &p2, bool e0, bool e1, bool e2, double z, bool mirrored);
	};


	class SpherePrimitive : public VboPrimitive
	{
		public:
			SpherePrimitive(float radius, bool center, Operation operation, unsigned int convexity, bool use_vaos);
			virtual ~SpherePrimitive() {};
			virtual void initialize();
		private:
			float radius_;
			bool center_;
	};

	class CylinderPrimitive : public VboPrimitive
	{
		public:
			CylinderPrimitive(float radius1, float radius2, float height, bool center, Operation operation, unsigned int convexity, bool use_vaos);
			virtual ~CylinderPrimitive() {};
			virtual void initialize();
		private:
			float radius1_;
			float radius2_;
			float height_;
			bool center_;
	};

	class CubePrimitive : public VboPrimitive
	{
		public:
			CubePrimitive(float x, float y, float z, bool center, Operation operation, unsigned int convexity, bool use_vaos);
			virtual ~CubePrimitive() {};
			virtual void initialize();
		private:
			float x_, y_, z_;
			bool center_;
	};
} // namespace OpenCSG
