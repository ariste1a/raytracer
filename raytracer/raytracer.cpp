//// [header]
//// A very basic raytracer example.
//// [/header]
//// [compile]
//// c++ -o raytracer -O3 -Wall raytracer.cpp
//// [/compile]
//// [ignore]
//// Copyright (C) 2012  www.scratchapixel.com
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//// [/ignore]
//#include <cstdlib>
//#include <cstdio>
//#include <cmath>
//#include <fstream>
//#include <vector>
//#include <iostream>
//#include <cassert>
//#include "tiny_obj_loader.h"
//#include <windows.h>
//#include <gdiplus.h>
//
//#define TINYOBJLOADER_IMPLEMENTATION 
//#if defined __linux__ || defined __APPLE__
//// "Compiled for Linux
//#else
//// Windows doesn't define these values by default, Linux does
//#define M_PI 3.141592653589793
////#define INFINITY 1e8
//#endif
//
//using namespace Gdiplus;
//#pragma comment (lib,"Gdiplus.lib")
//template<typename T>
//class Vec3
//{
//public:
//	T x, y, z;
//	Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
//	Vec3(T xx) : x(xx), y(xx), z(xx) {}
//	Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
//	Vec3& normalize()
//	{
//		T nor2 = length2();
//		if (nor2 > 0) {
//			T invNor = 1 / sqrt(nor2);
//			x *= invNor, y *= invNor, z *= invNor;
//		}
//		return *this;
//	}
//	Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
//	Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
//	T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
//	Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
//	Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
//	Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
//	Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
//	Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
//	T length2() const { return x * x + y * y + z * z; }
//	T length() const { return sqrt(length2()); }
//	friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
//	{
//		os << "[" << v.x << " " << v.y << " " << v.z << "]";
//		return os;
//	}
//};
//
//typedef Vec3<float> Vec3f;
//
//typedef struct {
//	//GLuint vb;  // vertex buffer
//	int numTriangles;
//	size_t material_id;
//} DrawObject;
//
//class Sphere
//{
//public:
//	Vec3f center;                           /// position of the sphere
//	float radius, radius2;                  /// sphere radius and radius^2
//	Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
//	float transparency, reflection;         /// surface transparency and reflectivity
//	Sphere(
//		const Vec3f &c,
//		const float &r,
//		const Vec3f &sc,
//		const float &refl = 0,
//		const float &transp = 0,
//		const Vec3f &ec = 0) :
//		center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
//		transparency(transp), reflection(refl)
//	{ /* empty */
//	}
//	//[comment]
//	// Compute a ray-sphere intersection using the geometric solution
//	//[/comment]
//	bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
//	{
//		Vec3f l = center - rayorig;
//		float tca = l.dot(raydir);
//		if (tca < 0) return false;
//		float d2 = l.dot(l) - tca * tca;
//		if (d2 > radius2) return false;
//		float thc = sqrt(radius2 - d2);
//		t0 = tca - thc;
//		t1 = tca + thc;
//
//		return true;
//	}
//};
//
////[comment]
//// This variable controls the maximum recursion depth
////[/comment]
//#define MAX_RAY_DEPTH 5
//
//float mix(const float &a, const float &b, const float &mix)
//{
//	return b * mix + a * (1 - mix);
//}
//
////[comment]
//// This is the main trace function. It takes a ray as argument (defined by its origin
//// and direction). We test if this ray intersects any of the geometry in the scene.
//// If the ray intersects an object, we compute the intersection point, the normal
//// at the intersection point, and shade this point using this information.
//// Shading depends on the surface property (is it transparent, reflective, diffuse).
//// The function returns a color for the ray. If the ray intersects an object that
//// is the color of the object at the intersection point, otherwise it returns
//// the background color.
////[/comment]
//Vec3f trace(
//	const Vec3f &rayorig,
//	const Vec3f &raydir,
//	const std::vector<Sphere> &spheres,
//	const int &depth)
//{
//	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
//	float tnear = INFINITY;
//	const Sphere* sphere = NULL;
//	// find intersection of this ray with the sphere in the scene
//	for (unsigned i = 0; i < spheres.size(); ++i) {
//		float t0 = INFINITY, t1 = INFINITY;
//		if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
//			if (t0 < 0) t0 = t1;
//			if (t0 < tnear) {
//				tnear = t0;
//				sphere = &spheres[i];
//			}
//		}
//	}
//	// if there's no intersection return black or background color
//	if (!sphere) return Vec3f(2);
//	Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
//	Vec3f phit = rayorig + raydir * tnear; // point of intersection
//	Vec3f nhit = phit - sphere->center; // normal at the intersection point
//	nhit.normalize(); // normalize normal direction
//					  // If the normal and the view direction are not opposite to each other
//					  // reverse the normal direction. That also means we are inside the sphere so set
//					  // the inside bool to true. Finally reverse the sign of IdotN which we want
//					  // positive.
//	float bias = 1e-4; // add some bias to the point from which we will be tracing
//	bool inside = false;
//	if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
//	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
//		float facingratio = -raydir.dot(nhit);
//		// change the mix value to tweak the effect
//		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
//		// compute reflection direction (not need to normalize because all vectors
//		// are already normalized)
//		Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
//		refldir.normalize();
//		Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
//		Vec3f refraction = 0;
//		// if the sphere is also transparent compute refraction ray (transmission)
//		if (sphere->transparency) {
//			float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
//			float cosi = -nhit.dot(raydir);
//			float k = 1 - eta * eta * (1 - cosi * cosi);
//			Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
//			refrdir.normalize();
//			refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
//		}
//		// the result is a mix of reflection and refraction (if the sphere is transparent)
//		surfaceColor = (
//			reflection * fresneleffect +
//			refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
//	}
//	else {
//		// it's a diffuse object, no need to raytrace any further
//		for (unsigned i = 0; i < spheres.size(); ++i) {
//			if (spheres[i].emissionColor.x > 0) {
//				// this is a light
//				Vec3f transmission = 1;
//				Vec3f lightDirection = spheres[i].center - phit;
//				lightDirection.normalize();
//				for (unsigned j = 0; j < spheres.size(); ++j) {
//					if (i != j) {
//						float t0, t1;
//						if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
//							transmission = 0;
//							break;
//						}
//					}
//				}
//				surfaceColor += sphere->surfaceColor * transmission *
//					std::fmax(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
//			}
//		}
//	}
//
//	return surfaceColor + sphere->emissionColor;
//}
//
////[comment]
//// Main rendering function. We compute a camera ray for each pixel of the image
//// trace it and return a color. If the ray hits a sphere, we return the color of the
//// sphere at the intersection point, else we return the background color.
////[/comment]
//void render(const std::vector<Sphere> &spheres)
//{
//	unsigned width = 640, height = 480;
//	Vec3f *image = new Vec3f[width * height], *pixel = image;
//	float invWidth = 1 / float(width), invHeight = 1 / float(height);
//	float fov = 30, aspectratio = width / float(height);
//	float angle = tan(M_PI * 0.5 * fov / 180.);
//	// Trace rays
//	for (unsigned y = 0; y < height; ++y) {
//		for (unsigned x = 0; x < width; ++x, ++pixel) {
//			float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
//			float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
//			Vec3f raydir(xx, yy, -1);
//			raydir.normalize();
//			*pixel = trace(Vec3f(0), raydir, spheres, 0);
//		}
//	}
//	// Save result to a PPM image (keep these flags if you compile under Windows)
//	std::ofstream ofs("./untitled.bmp", std::ios::out | std::ios::binary);
//	//ofs << "P6\n" << width << " " << height << "\n255\n";
//	for (unsigned i = 0; i < width * height; ++i) {
//		ofs << (unsigned char)(std::fmin(float(1), image[i].x) * 255) <<
//			(unsigned char)(std::fmin(float(1), image[i].y) * 255) <<
//			(unsigned char)(std::fmin(float(1), image[i].z) * 255);
//	}
//	
//	for (unsigned i = 0; i < width * height; ++i) {
//		ofs << (unsigned char)(std::fmin(float(1), image[i].x) * 255) <<
//			(unsigned char)(std::fmin(float(1), image[i].y) * 255) <<
//			(unsigned char)(std::fmin(float(1), image[i].z) * 255);
//	}	
//	ofs.close();
//	delete[] image;
//}
//
//void loadObj()
//{
//	std::string inputfile = "CornellBox-Original.obj";
//	tinyobj::attrib_t attrib;
//	std::vector<tinyobj::shape_t> shapes;
//	std::vector<tinyobj::material_t> materials;
//
//	std::string err;
//	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());
//
//	if (!err.empty()) { // `err` may contain warning message.
//		std::cerr << err << std::endl;
//	}
//
//	if (!ret) {
//		exit(1);
//	}
//
//	// Loop over shapes
//	for (size_t s = 0; s < shapes.size(); s++) {
//		// Loop over faces(polygon)
//		size_t index_offset = 0;
//		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
//			int fv = shapes[s].mesh.num_face_vertices[f];
//
//			// Loop over vertices in the face.
//			for (size_t v = 0; v < fv; v++) {
//				// access to vertex
//				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
//				tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
//				tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
//				tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
//				tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
//				tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
//				tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
//				tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
//				tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];
//			}
//			index_offset += fv;
//
//			// per-face material
//			shapes[s].mesh.material_ids[f];
//		}
//	}
//}
//
//static std::string GetBaseDir(const std::string &filepath) {
//	if (filepath.find_last_of("/\\") != std::string::npos)
//		return filepath.substr(0, filepath.find_last_of("/\\"));
//	return "";
//}
//
//static bool FileExists(const std::string &abs_filename) {
//	bool ret;
//	FILE *fp = fopen(abs_filename.c_str(), "rb");
//	if (fp) {
//		ret = true;
//		fclose(fp);
//	}
//	else {
//		ret = false;
//	}
//
//	return ret;
//}
//
//static void CalcNormal(float N[3], float v0[3], float v1[3], float v2[3]) {
//	float v10[3];
//	v10[0] = v1[0] - v0[0];
//	v10[1] = v1[1] - v0[1];
//	v10[2] = v1[2] - v0[2];
//
//	float v20[3];
//	v20[0] = v2[0] - v0[0];
//	v20[1] = v2[1] - v0[1];
//	v20[2] = v2[2] - v0[2];
//
//	N[0] = v20[1] * v10[2] - v20[2] * v10[1];
//	N[1] = v20[2] * v10[0] - v20[0] * v10[2];
//	N[2] = v20[0] * v10[1] - v20[1] * v10[0];
//
//	float len2 = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
//	if (len2 > 0.0f) {
//		float len = sqrtf(len2);
//
//		N[0] /= len;
//		N[1] /= len;
//	}
//}
//#undef max
//#undef min
//static bool LoadObjAndConvert(float bmin[3], float bmax[3],
//	std::vector<DrawObject>* drawObjects,
//	std::vector<tinyobj::material_t>& materials,
//	std::map<std::string, Bitmap*>& textures,
//	const char* filename) {
//	tinyobj::attrib_t attrib;
//	std::vector<tinyobj::shape_t> shapes;
//
//	//timerutil tm;
//
//	//tm.start();
//
//	std::string base_dir = GetBaseDir(filename);
//#ifdef _WIN32
//	base_dir += "\\";
//#else
//	base_dir += "/";
//#endif
//
//	std::string err;
//	bool ret =
//		tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename, base_dir.c_str());
//	if (!err.empty()) {
//		std::cerr << err << std::endl;
//	}
//
//	// tm.end();
//
//	if (!ret) {
//		std::cerr << "Failed to load " << filename << std::endl;
//		return false;
//	}
//
//	//printf("Parsing time: %d [ms]\n", (int)tm.msec());
//
//	printf("# of vertices  = %d\n", (int)(attrib.vertices.size()) / 3);
//	printf("# of normals   = %d\n", (int)(attrib.normals.size()) / 3);
//	printf("# of texcoords = %d\n", (int)(attrib.texcoords.size()) / 2);
//	printf("# of materials = %d\n", (int)materials.size());
//	printf("# of shapes    = %d\n", (int)shapes.size());
//
//	// Append `default` material
//	materials.push_back(tinyobj::material_t());
//
//	// Load diffuse textures
//	{
//		for (size_t m = 0; m < materials.size(); m++) {
//			tinyobj::material_t* mp = &materials[m];
//
//			if (mp->diffuse_texname.length() > 0) {
//				// Only load the texture if it is not already loaded
//				if (textures.find(mp->diffuse_texname) == textures.end()) {
//					//GLuint texture_id;
//					int w, h;
//					int comp;
//
//					std::string texture_filename = mp->diffuse_texname;
//					if (!FileExists(texture_filename)) {
//						// Append base dir.
//						texture_filename = base_dir + mp->diffuse_texname;
//						if (!FileExists(texture_filename)) {
//							std::cerr << "Unable to find file: " << mp->diffuse_texname << std::endl;
//							exit(1);
//						}
//					}
//					std::wstring name(L"texture_filename");
//					Bitmap* bitmap = Bitmap::FromFile(name.c_str());
//					/*unsigned char* image = stbi_load(texture_filename.c_str(), &w, &h, &comp, STBI_default);
//					if (!image) {
//						std::cerr << "Unable to load texture: " << texture_filename << std::endl;
//						exit(1);
//					}
//					glGenTextures(1, &texture_id);
//					glBindTexture(GL_TEXTURE_2D, texture_id);
//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//					if (comp == 3) {
//						glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
//					}
//					else if (comp == 4) {
//						glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
//					}
//					glBindTexture(GL_TEXTURE_2D, 0);
//					stbi_image_free(image);*/
//					textures.insert(std::make_pair(mp->diffuse_texname, bitmap));
//				}
//			}
//		}
//	}
//
//	bmin[0] = bmin[1] = bmin[2] = std::numeric_limits<float>::max(); 
//	bmax[0] = bmax[1] = bmax[2] = -std::numeric_limits<float>::max();
//
//	{
//		for (size_t s = 0; s < shapes.size(); s++) {
//			DrawObject o;
//			std::vector<float> vb;  // pos(3float), normal(3float), color(3float)
//			for (size_t f = 0; f < shapes[s].mesh.indices.size() / 3; f++) {
//				tinyobj::index_t idx0 = shapes[s].mesh.indices[3 * f + 0];
//				tinyobj::index_t idx1 = shapes[s].mesh.indices[3 * f + 1];
//				tinyobj::index_t idx2 = shapes[s].mesh.indices[3 * f + 2];
//
//				int current_material_id = shapes[s].mesh.material_ids[f];
//
//				if ((current_material_id < 0) || (current_material_id >= static_cast<int>(materials.size()))) {
//					// Invaid material ID. Use default material.
//					current_material_id = materials.size() - 1; // Default material is added to the last item in `materials`.
//				}
//				//if (current_material_id >= materials.size()) {
//				//    std::cerr << "Invalid material index: " << current_material_id << std::endl;
//				//}
//				//
//				float diffuse[3];
//				for (size_t i = 0; i < 3; i++) {
//					diffuse[i] = materials[current_material_id].diffuse[i];
//				}
//				float tc[3][2];
//				if (attrib.texcoords.size() > 0) {
//					assert(attrib.texcoords.size() > 2 * idx0.texcoord_index + 1);
//					assert(attrib.texcoords.size() > 2 * idx1.texcoord_index + 1);
//					assert(attrib.texcoords.size() > 2 * idx2.texcoord_index + 1);
//					tc[0][0] = attrib.texcoords[2 * idx0.texcoord_index];
//					tc[0][1] = 1.0f - attrib.texcoords[2 * idx0.texcoord_index + 1];
//					tc[1][0] = attrib.texcoords[2 * idx1.texcoord_index];
//					tc[1][1] = 1.0f - attrib.texcoords[2 * idx1.texcoord_index + 1];
//					tc[2][0] = attrib.texcoords[2 * idx2.texcoord_index];
//					tc[2][1] = 1.0f - attrib.texcoords[2 * idx2.texcoord_index + 1];
//				}
//				else {
//					tc[0][0] = 0.0f;
//					tc[0][1] = 0.0f;
//					tc[1][0] = 0.0f;
//					tc[1][1] = 0.0f;
//					tc[2][0] = 0.0f;
//					tc[2][1] = 0.0f;
//				}
//
//				float v[3][3];
//				for (int k = 0; k < 3; k++) {
//					int f0 = idx0.vertex_index;
//					int f1 = idx1.vertex_index;
//					int f2 = idx2.vertex_index;
//					assert(f0 >= 0);
//					assert(f1 >= 0);
//					assert(f2 >= 0);
//
//					v[0][k] = attrib.vertices[3 * f0 + k];
//					v[1][k] = attrib.vertices[3 * f1 + k];
//					v[2][k] = attrib.vertices[3 * f2 + k];
//					bmin[k] = std::fmin(v[0][k], bmin[k]);
//					bmin[k] = std::fmin(v[1][k], bmin[k]);
//					bmin[k] = std::fmin(v[2][k], bmin[k]);
//					bmax[k] = std::fmax(v[0][k], bmax[k]);
//					bmax[k] = std::fmax(v[1][k], bmax[k]);
//					bmax[k] = std::fmax(v[2][k], bmax[k]);
//				}
//
//				float n[3][3];
//				if (attrib.normals.size() > 0) {
//					int f0 = idx0.normal_index;
//					int f1 = idx1.normal_index;
//					int f2 = idx2.normal_index;
//					assert(f0 >= 0);
//					assert(f1 >= 0);
//					assert(f2 >= 0);
//					for (int k = 0; k < 3; k++) {
//						n[0][k] = attrib.normals[3 * f0 + k];
//						n[1][k] = attrib.normals[3 * f1 + k];
//						n[2][k] = attrib.normals[3 * f2 + k];
//					}
//				}
//				else {
//					// compute geometric normal
//					CalcNormal(n[0], v[0], v[1], v[2]);
//					n[1][0] = n[0][0];
//					n[1][1] = n[0][1];
//					n[1][2] = n[0][2];
//					n[2][0] = n[0][0];
//					n[2][1] = n[0][1];
//					n[2][2] = n[0][2];
//				}
//
//				for (int k = 0; k < 3; k++) {
//					vb.push_back(v[k][0]);
//					vb.push_back(v[k][1]);
//					vb.push_back(v[k][2]);
//					vb.push_back(n[k][0]);
//					vb.push_back(n[k][1]);
//					vb.push_back(n[k][2]);
//					// Combine normal and diffuse to get color.
//					float normal_factor = 0.2;
//					float diffuse_factor = 1 - normal_factor;
//					float c[3] = {
//						n[k][0] * normal_factor + diffuse[0] * diffuse_factor,
//						n[k][1] * normal_factor + diffuse[1] * diffuse_factor,
//						n[k][2] * normal_factor + diffuse[2] * diffuse_factor
//					};
//					float len2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
//					if (len2 > 0.0f) {
//						float len = sqrtf(len2);
//
//						c[0] /= len;
//						c[1] /= len;
//						c[2] /= len;
//					}
//					vb.push_back(c[0] * 0.5 + 0.5);
//					vb.push_back(c[1] * 0.5 + 0.5);
//					vb.push_back(c[2] * 0.5 + 0.5);
//
//					vb.push_back(tc[k][0]);
//					vb.push_back(tc[k][1]);
//				}
//			}
//
//			//o.vb = 0;
//			o.numTriangles = 0;
//
//			// OpenGL viewer does not support texturing with per-face material.
//			if (shapes[s].mesh.material_ids.size() > 0 && shapes[s].mesh.material_ids.size() > s) {
//				// Base case
//				o.material_id = shapes[s].mesh.material_ids[s];
//			}
//			else {
//				o.material_id = materials.size() - 1; // = ID for default material.
//			}
//
//			if (vb.size() > 0) {
//				/*glGenBuffers(1, &o.vb);
//				glBindBuffer(GL_ARRAY_BUFFER, o.vb);
//				glBufferData(GL_ARRAY_BUFFER, vb.size() * sizeof(float), &vb.at(0),
//					GL_STATIC_DRAW);*/
//				o.numTriangles = vb.size() / (3 + 3 + 3 + 2) / 3; // 3:vtx, 3:normal, 3:col, 2:texcoord
//
//				printf("shape[%d] # of triangles = %d\n", static_cast<int>(s),
//					o.numTriangles);
//			}
//
//			drawObjects->push_back(o);
//		}
//	}
//
//	printf("bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
//	printf("bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
//
//	return true;
//}
//
////[comment]
//// In the main function, we will create the scene which is composed of 5 spheres
//// and 1 light (which is also a sphere). Then, once the scene description is complete
//// we render that scene, by calling the render() function.
////[/comment]
//
//std::vector<DrawObject> gDrawObjects;
//
////int main(int argc, char **argv)
////{
////	//HWND                hWnd;
////	//MSG                 msg;
////	//WNDCLASS            wndClass;
////	//GdiplusStartupInput gdiplusStartupInput;
////	//ULONG_PTR           gdiplusToken;
////
////	//// Initialize GDI+.
////	//GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
////
////
////	srand(13);
////	std::vector<Sphere> spheres;
////	//loadObj(); 
////	float bmin[3], bmax[3]; 
////	std::vector<tinyobj::material_t> materials;
////	std::map<std::string, Bitmap*> textures;
////	LoadObjAndConvert(bmin, bmax, &gDrawObjects, materials, textures, "./CornellBox-Original.obj");
////	// position, radius, surface color, reflectivity, transparency, emission color
////	spheres.push_back(Sphere(Vec3f(0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
////	spheres.push_back(Sphere(Vec3f(0.0, 0, -20), 4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
////	spheres.push_back(Sphere(Vec3f(5.0, -1, -15), 2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
////	spheres.push_back(Sphere(Vec3f(5.0, 0, -25), 3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
////	spheres.push_back(Sphere(Vec3f(-5.5, 0, -15), 3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
////	// light
////	spheres.push_back(Sphere(Vec3f(0.0, 20, -30), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
////	render(spheres);
////
////	//GdiplusShutdown(gdiplusToken);	
////	return 0;
////}