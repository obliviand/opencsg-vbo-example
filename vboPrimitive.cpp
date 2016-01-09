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
#include <GL/glew.h>
#include "vboPrimitive.h"

namespace OpenCSG {

	VboPrimitive::VboPrimitive(Operation operation, unsigned int convexity, bool use_vaos)
		: Primitive(operation, convexity), vao_id_(0), vbo_id_(0), triangle_count_(0), use_vaos_(use_vaos)
	{
		if (use_vaos) {
			glGenVertexArrays(1, &vao_id_);
		}
		checkGLErrors("render_surface glGenVertexArrays");
		m_ = Transform3d::Identity();
	}

	VboPrimitive::~VboPrimitive()
	{
		if (use_vaos_) {
			glDeleteVertexArrays(1, &vao_id_);
		} else {
			glDeleteBuffers(1, &vbo_id_);
		}
	}

	GLuint VboPrimitive::getVaoId() const {
		return vao_id_;
	}

	void VboPrimitive::rotate(float x, float y, float z)
	{
		Eigen::AngleAxisd rotx(0, Vector3d::UnitX());
		Eigen::AngleAxisd roty(0, Vector3d::UnitY());
		Eigen::AngleAxisd rotz(0, Vector3d::UnitZ());
		rotx = Eigen::AngleAxisd(x*M_PI/180, Vector3d::UnitX());
		roty = Eigen::AngleAxisd(y*M_PI/180, Vector3d::UnitY());
		rotz = Eigen::AngleAxisd(z*M_PI/180, Vector3d::UnitZ());
		m_.rotate(rotz * roty * rotx);
	}

	void VboPrimitive::translate(float x, float y, float z)
	{
		Vector3d translatevec(x,y,z);
		m_.translate(translatevec);
	}

	void VboPrimitive::scale(float x, float y, float z)
	{
		Vector3d scalevec(x,y,z);
		m_.scale(scalevec);
	}

	void VboPrimitive::render()
	{
		GLint current_program;
		GLint current_vao;
		GLint current_vbo;

		glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);
		glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &current_vao);
		glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &current_vbo);

		glPushMatrix();
			glMultMatrixd(m_.data());

			GLboolean origVertexArrayState = glIsEnabled(GL_VERTEX_ARRAY);
			GLboolean origNormalArrayState = glIsEnabled(GL_NORMAL_ARRAY);
			if (!origVertexArrayState) {
				glEnableClientState(GL_VERTEX_ARRAY);
			}
			if (!origNormalArrayState) {
				glEnableClientState(GL_NORMAL_ARRAY);
			}

			if (use_vaos_) {
				glBindVertexArray(vao_id_);
			} else {
				glBindBuffer(GL_ARRAY_BUFFER, vbo_id_);
				glVertexPointer(3, GL_FLOAT, sizeof(Vertex3d), (GLvoid *)(0));
				glNormalPointer(GL_FLOAT, sizeof(Vertex3d), (GLvoid *)(sizeof(float)*3));
			}

				glDrawArrays(GL_TRIANGLES, 0, triangle_count_);
			if (use_vaos_) {
				glBindVertexArray(current_vao);
			} else {
				glBindBuffer(GL_ARRAY_BUFFER, current_vbo);
			}

			if (!origVertexArrayState) {
				glDisableClientState(GL_VERTEX_ARRAY);
			}
			if (!origNormalArrayState) {
				glDisableClientState(GL_NORMAL_ARRAY);
			}
		glPopMatrix();
	}

	void VboPrimitive::draw_triangle(Vertex3d *vertices,
					 const Vertex3d &p0, const Vertex3d &p1, const Vertex3d &p2,
					 double e0f, double e1f, double e2f,
					 double nx, double ny, double nz, double nl,
					 double z, bool mirror)
	{
		// vertex
		vertices[0].x = p0.x;
		vertices[0].y = p0.y;
		vertices[0].z = p0.z + z;
		// normal
		vertices[0].nx = nx/nl;
		vertices[0].ny = ny/nl;
		vertices[0].nz = nz/nl;
		// color
		vertices[0].r = 0.0;
		vertices[0].g = 0.0;
		vertices[0].b = 1.0;

		if (!mirror) {
			//vertex
			vertices[1].x = p1.x;
			vertices[1].y = p1.y;
			vertices[1].z = p1.z + z;
			// normal
			vertices[1].nx = nx/nl;
			vertices[1].ny = ny/nl;
			vertices[1].nz = nz/nl;
			// color
			vertices[1].r = 0.0;
			vertices[1].g = 0.0;
			vertices[1].b = 1.0;

			// vertex
			vertices[2].x = p2.x;
			vertices[2].y = p2.y;
			vertices[2].z = p2.z + z;
			// normal
			vertices[2].nx = nx/nl;
			vertices[2].ny = ny/nl;
			vertices[2].nz = nz/nl;
			// color
			vertices[2].r = 0.0;
			vertices[2].g = 0.0;
			vertices[2].b = 1.0;
		}


		if (mirror) {
			// vertex
			vertices[1].x = p2.x;
			vertices[1].y = p2.y;
			vertices[1].z = p2.z + z;
			// normal
			vertices[1].nx = nx/nl;
			vertices[1].ny = ny/nl;
			vertices[1].nz = nz/nl;
			// color
			vertices[1].r = 0.0;
			vertices[1].g = 0.0;
			vertices[1].b = 1.0;

			// vertex
			vertices[2].x = p1.x;
			vertices[2].y = p1.y;
			vertices[2].z = p1.z + z;
			// normal
			vertices[2].nx = nx/nl;
			vertices[2].ny = ny/nl;
			vertices[2].nz = nz/nl;
			// color
			vertices[2].r = 0.0;
			vertices[2].g = 0.0;
			vertices[2].b = 1.0;
		}
	}

	void VboPrimitive::gl_draw_triangle(Vertex3d *vertices, const Vertex3d &p0, const Vertex3d &p1, const Vertex3d &p2, bool e0, bool e1, bool e2, double z, bool mirrored)
	{
		double ax = p1.x - p0.x, bx = p1.x - p2.x;
		double ay = p1.y - p0.y, by = p1.y - p2.y;
		double az = p1.z - p0.z, bz = p1.z - p2.z;
		double nx = ay*bz - az*by;
		double ny = az*bx - ax*bz;
		double nz = ax*by - ay*bx;
		double nl = sqrt(nx*nx + ny*ny + nz*nz);

		double e0f = 0.0;
		double e1f = 0.0;
		double e2f = 0.0;

		draw_triangle(vertices, p0, p1, p2, e0f, e1f, e2f, nx, ny, nz, nl, z, mirrored);
	}

	void VboPrimitive::render_surface(const Polygons &polygons)
	{
		bool mirrored = m_.matrix().determinant() < 0;
		size_t buffer_offset = 0;
		size_t buffer_size = 0;
		VboPrimitive::Vertex3d *vertices;
		// pre-calc size
		for (size_t i = 0; i < polygons.size(); i++) {
			const Polygon *poly = &polygons[i];
			if (poly->size() == 3 || poly->size() > 4) {
				buffer_size += poly->size()*sizeof(Vertex3d)*3; // One triangle per polygon
			} else if (poly->size() == 4) {
				buffer_size += 2*sizeof(Vertex3d)*3; // Two triangles per polygon
			}
		}
		vertices = new Vertex3d[buffer_size];

		for (size_t i = 0; i < polygons.size(); i++) {
			const Polygon *poly = &polygons[i];

			if (poly->size() == 3) {
				gl_draw_triangle(&vertices[buffer_offset], poly->at(0), poly->at(1), poly->at(2), true, true, true, 0, mirrored);
				buffer_offset+=3;
			} else if (poly->size() == 4) {
				gl_draw_triangle(&vertices[buffer_offset], poly->at(0), poly->at(1), poly->at(3), true, false, true, 0, mirrored);
				buffer_offset+=3;
				gl_draw_triangle(&vertices[buffer_offset], poly->at(2), poly->at(3), poly->at(1), true, false, true, 0, mirrored);
				buffer_offset+=3;
			} else {
				Vertex3d center;
				center.x = 0.0; center.y = 0.0; center.z = 0.0;
				for (size_t j = 0; j < poly->size(); j++) {
					center.x += poly->at(j).x;
					center.y += poly->at(j).y;
					center.z += poly->at(j).z;
				}
				center.x /= poly->size();
				center.y /= poly->size();
				center.z /= poly->size();
				for (size_t j = 1; j <= poly->size(); j++) {
					gl_draw_triangle(&vertices[buffer_offset], center, poly->at(j - 1), poly->at(j % poly->size()), false, true, false, 0, mirrored);
					buffer_offset+=3;
				}
			}
		}

		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindVertexArray(vao_id_);
			checkGLErrors("render_surface glBindVertexArray");
			glBindBuffer(GL_ARRAY_BUFFER, vbo);
				glBufferData(GL_ARRAY_BUFFER, buffer_size, vertices, GL_STATIC_DRAW);
				glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex3d), (void *)(0));
				glEnableVertexAttribArray(0);
				glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex3d), (void *)(sizeof(float)*3));
				glEnableVertexAttribArray(2);
		glBindVertexArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		if (use_vaos_) {
			glDeleteBuffers(1, &vbo);
		} else {
			vbo_id_ = vbo;
		}

		triangle_count_ = buffer_size / sizeof(Vertex3d);
		delete [] vertices;
	}

	SpherePrimitive::SpherePrimitive(float radius, bool center, Operation operation, unsigned int convexity, bool use_vaos)
		: VboPrimitive(operation, convexity, use_vaos), radius_(radius), center_(center)
	{
	}

	void SpherePrimitive::initialize()
	{
		Polygons poly;
		if (radius_ > 0 && !isinf(radius_)) {
			struct ring_s {
				Point2d *points;
				double z;
			};

			int fragments = get_fragments_from_r(radius_);
			int rings = (fragments+1)/2;
// Uncomment the following three lines to enable experimental sphere tesselation
//		if (rings % 2 == 0) rings++; // To ensure that the middle ring is at phi == 0 degrees

			ring_s *ring = new ring_s[rings];

//		double offset = 0.5 * ((fragments / 2) % 2);
			for (int i = 0; i < rings; i++) {
//			double phi = (M_PI * (i + offset)) / (fragments/2);
				double phi = (M_PI * (i + 0.5)) / rings;
				double r = radius_ * sin(phi);
				ring[i].z = radius_ * cos(phi);
				ring[i].points = new Point2d[fragments];
				generate_circle(ring[i].points, r, fragments);
			}

			poly.push_back(Polygon());
			for (int i = 0; i < fragments; i++) {
				Vertex3d vertex;
				vertex.x = ring[0].points[i].x;
				vertex.y = ring[0].points[i].y;
				vertex.z = ring[0].z;
				poly.back().push_back(vertex);
			}

			for (int i = 0; i < rings-1; i++) {
				ring_s *r1 = &ring[i];
				ring_s *r2 = &ring[i+1];
				int r1i = 0, r2i = 0;
				Vertex3d vertex;
				while (r1i < fragments || r2i < fragments)
				{
					if (r1i >= fragments)
						goto sphere_next_r2;
					if (r2i >= fragments)
						goto sphere_next_r1;
					if ((double)r1i / fragments < (double)r2i / fragments)
					{
						sphere_next_r1:
							poly.push_back(Polygon());
							int r1j = (r1i+1) % fragments;
							vertex.x = r1->points[r1i].x;
							vertex.y = r1->points[r1i].y;
							vertex.z = r1->z;
							poly.back().insert(poly.back().begin(),vertex);
							vertex.x = r1->points[r1j].x;
							vertex.y = r1->points[r1j].y;
							vertex.z = r1->z;
							poly.back().insert(poly.back().begin(),vertex);
							vertex.x = r2->points[r2i % fragments].x;
							vertex.y = r2->points[r2i % fragments].y;
							vertex.z = r2->z;
							poly.back().insert(poly.back().begin(),vertex);
							r1i++;
					} else {
						sphere_next_r2:
							poly.push_back(Polygon());
							int r2j = (r2i+1) % fragments;
							vertex.x = r2->points[r2i].x;
							vertex.y = r2->points[r2i].y;
							vertex.z = r2->z;
							poly.back().push_back(vertex);
							vertex.x = r2->points[r2j].x;
							vertex.y = r2->points[r2j].y;
							vertex.z = r2->z;
							poly.back().push_back(vertex);
							vertex.x = r1->points[r1i % fragments].x;
							vertex.y = r1->points[r1i % fragments].y;
							vertex.z = r1->z;
							poly.back().push_back(vertex);
							r2i++;
					}
				}
			}

			poly.push_back(Polygon());
			for (int i = 0; i < fragments; i++) {
				Vertex3d vertex;
				vertex.x = ring[rings-1].points[i].x;
				vertex.y = ring[rings-1].points[i].y;
				vertex.z = ring[rings-1].z;
				poly.back().insert(poly.back().begin(),vertex);
			}

			for (int i = 0; i < rings; i++) {
				delete[] ring[i].points;
			}
			delete[] ring;
		}
		render_surface(poly);
	}


	CylinderPrimitive::CylinderPrimitive(float radius1, float radius2, float height, bool center, Operation operation, unsigned int convexity, bool use_vaos)
		: VboPrimitive(operation, convexity, use_vaos), radius1_(radius1), radius2_(radius2), height_(height), center_(center)
	{
	}

	void CylinderPrimitive::initialize()
	{
		Polygons poly;
		if (height_ > 0 && !isinf(height_) &&
		    radius1_ >=0 && radius2_ >= 0 && (radius1_ > 0 || radius2_ > 0) &&
		    !isinf(radius1_) && !isinf(radius2_)) {
			int fragments = get_fragments_from_r(fmax(radius1_, radius2_));

			double z1, z2;
			if (center_) {
				z1 = -height_/2;
				z2 = +height_/2;
			} else {
				z1 = 0;
				z2 = height_;
			}

			Point2d *circle1 = new Point2d[fragments];
			Point2d *circle2 = new Point2d[fragments];

			generate_circle(circle1, radius1_, fragments);
			generate_circle(circle2, radius2_, fragments);

			for (int i=0; i<fragments; i++) {
				int j = (i+1) % fragments;
				if (radius1_ == radius2_) {
					Vertex3d vertex;
					poly.push_back(Polygon());
					vertex.x = circle1[i].x; vertex.y = circle1[i].y; vertex.z = z1;
					poly.back().insert(poly.back().begin(),vertex);
					vertex.x = circle2[i].x; vertex.y = circle2[i].y; vertex.z = z2;
					poly.back().insert(poly.back().begin(),vertex);
					vertex.x = circle2[j].x; vertex.y = circle2[j].y; vertex.z = z2;
					poly.back().insert(poly.back().begin(),vertex);
					vertex.x = circle1[j].x; vertex.y = circle1[j].y; vertex.z = z1;
					poly.back().insert(poly.back().begin(),vertex);
				} else {
					if (radius1_ > 0) {
						Vertex3d vertex;
						poly.push_back(Polygon());
						vertex.x = circle1[i].x; vertex.y = circle1[i].y; vertex.z = z1;
						poly.back().insert(poly.back().begin(),vertex);
						vertex.x = circle2[i].x; vertex.y = circle2[i].y; vertex.z = z2;
						poly.back().insert(poly.back().begin(),vertex);
						vertex.x = circle1[j].x; vertex.y = circle1[j].y; vertex.z = z1;
						poly.back().insert(poly.back().begin(),vertex);
					}
					if (radius2_ > 0) {
						Vertex3d vertex;
						poly.push_back(Polygon());
						vertex.x = circle2[i].x; vertex.y = circle2[i].y; vertex.z = z2;
						poly.back().insert(poly.back().begin(),vertex);
						vertex.x = circle2[j].x; vertex.y = circle2[j].y; vertex.z = z2;
						poly.back().insert(poly.back().begin(),vertex);
						vertex.x = circle1[j].x; vertex.y = circle1[j].y; vertex.z = z1;
						poly.back().insert(poly.back().begin(),vertex);
					}
				}
			}

			if (radius1_ > 0) {
				Vertex3d vertex;
				poly.push_back(Polygon());
				for (int i=0; i<fragments; i++) {
					vertex.x = circle1[i].x; vertex.y = circle1[i].y; vertex.z = z1;
					poly.back().insert(poly.back().begin(),vertex);
				}
			}

			if (radius2_ > 0) {
				Vertex3d vertex;
				poly.push_back(Polygon());
				for (int i=0; i<fragments; i++) {
					vertex.x = circle2[i].x; vertex.y = circle2[i].y; vertex.z = z2;
					poly.back().push_back(vertex);
				}
			}

			delete[] circle1;
			delete[] circle2;
		}
		render_surface(poly);
	}

	CubePrimitive::CubePrimitive(float x, float y, float z, bool center, Operation operation, unsigned int convexity, bool use_vaos)
		: VboPrimitive(operation, convexity, use_vaos), x_(x), y_(y), z_(z), center_(center)
	{
	}

	void CubePrimitive::initialize()
	{
		Polygons poly;
		if (x_ > 0 && y_ > 0 && z_ > 0 &&
		    !isinf(x_) && !isinf(y_) && !isinf(z_)) {
			double x1, x2, y1, y2, z1, z2;
			if (center_) {
				x1 = -x_/2;
				x2 = +x_/2;
				y1 = -y_/2;
				y2 = +y_/2;
				z1 = -z_/2;
				z2 = +z_/2;
			} else {
				x1 = y1 = z1 = 0;
				x2 = x_;
				y2 = y_;
				z2 = z_;
			}

			Vertex3d vertex;
			poly.push_back(Polygon()); // top
			vertex.x = x1; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);

			poly.push_back(Polygon()); // bottom
			vertex.x = x1; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);

			poly.push_back(Polygon()); // side1
			vertex.x = x1; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);

			poly.push_back(Polygon()); // side2
			vertex.x = x2; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);

			poly.push_back(Polygon()); // side3
			vertex.x = x2; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x2; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);

			poly.push_back(Polygon()); // side4
			vertex.x = x1; vertex.y = y2, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y1, vertex.z = z1;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y1, vertex.z = z2;
			poly.back().push_back(vertex);
			vertex.x = x1; vertex.y = y2, vertex.z = z2;
			poly.back().push_back(vertex);
		}
		render_surface(poly);
	}
}
