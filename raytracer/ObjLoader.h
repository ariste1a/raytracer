#pragma once

#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <vector> 
#include <windows.h>
#include <gdiplus.h>
#include <cassert>
#include <iostream>
#include "tiny_obj_loader.h"
#include "Vectors.h"
#pragma comment (lib,"Gdiplus.lib")
#define TINYOBJLOADER_IMPLEMENTATION 
using namespace Gdiplus;
class ObjLoader
{
public:
	ObjLoader();
	typedef struct {
		//GLuint vb;  // vertex buffer
		std::vector<uint32_t> indices;
		std::vector<Vec3f> vertices; 
		std::vector<Vec2f> texcoords;
		int numTriangles;
		size_t material_id;
	} DrawObject;

	static std::string GetBaseDir(const std::string &filepath) {
		if (filepath.find_last_of("/\\") != std::string::npos)
			return filepath.substr(0, filepath.find_last_of("/\\"));
		return "";
	}

	static bool FileExists(const std::string &abs_filename) {
		bool ret;
		FILE *fp = fopen(abs_filename.c_str(), "rb");
		if (fp) {
			ret = true;
			fclose(fp);
		}
		else {
			ret = false;
		}

		return ret;
	}

	static void CalcNormal(float N[3], float v0[3], float v1[3], float v2[3]) {
		float v10[3];
		v10[0] = v1[0] - v0[0];
		v10[1] = v1[1] - v0[1];
		v10[2] = v1[2] - v0[2];

		float v20[3];
		v20[0] = v2[0] - v0[0];
		v20[1] = v2[1] - v0[1];
		v20[2] = v2[2] - v0[2];

		N[0] = v20[1] * v10[2] - v20[2] * v10[1];
		N[1] = v20[2] * v10[0] - v20[0] * v10[2];
		N[2] = v20[0] * v10[1] - v20[1] * v10[0];

		float len2 = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
		if (len2 > 0.0f) {
			float len = sqrtf(len2);

			N[0] /= len;
			N[1] /= len;
		}
	}
#undef max
#undef min

	static bool ObjLoader::LoadObjAndConvert(float bmin[3], float bmax[3],
		std::vector<DrawObject>* drawObjects,
		std::vector<tinyobj::material_t>& materials,
		std::map<std::string, Bitmap*>& textures,
		const char* filename) {
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;

		//timerutil tm;

		//tm.start();

		std::string base_dir = GetBaseDir(filename);
#ifdef _WIN32
		base_dir += "\\";
#else
		base_dir += "/";
#endif

		std::string err;
		bool ret =
			tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename, base_dir.c_str());
		if (!err.empty()) {
			std::cerr << err << std::endl;
		}

		// tm.end();

		if (!ret) {
			std::cerr << "Failed to load " << filename << std::endl;
			return false;
		}

		//printf("Parsing time: %d [ms]\n", (int)tm.msec());

		printf("# of vertices  = %d\n", (int)(attrib.vertices.size()) / 3);
		printf("# of normals   = %d\n", (int)(attrib.normals.size()) / 3);
		printf("# of texcoords = %d\n", (int)(attrib.texcoords.size()) / 2);
		printf("# of materials = %d\n", (int)materials.size());
		printf("# of shapes    = %d\n", (int)shapes.size());

		// Append `default` material
		materials.push_back(tinyobj::material_t());

		// Load diffuse textures
		{
			for (size_t m = 0; m < materials.size(); m++) {
				tinyobj::material_t* mp = &materials[m];

				if (mp->diffuse_texname.length() > 0) {
					// Only load the texture if it is not already loaded
					if (textures.find(mp->diffuse_texname) == textures.end()) {
						//GLuint texture_id;
						int w, h;
						int comp;

						std::string texture_filename = mp->diffuse_texname;
						if (!FileExists(texture_filename)) {
							// Append base dir.
							texture_filename = base_dir + mp->diffuse_texname;
							if (!FileExists(texture_filename)) {
								std::cerr << "Unable to find file: " << mp->diffuse_texname << std::endl;
								exit(1);
							}
						}
						std::wstring name(L"texture_filename");
						Bitmap* bitmap = Bitmap::FromFile(name.c_str());
						/*unsigned char* image = stbi_load(texture_filename.c_str(), &w, &h, &comp, STBI_default);
						if (!image) {
						std::cerr << "Unable to load texture: " << texture_filename << std::endl;
						exit(1);
						}
						glGenTextures(1, &texture_id);
						glBindTexture(GL_TEXTURE_2D, texture_id);
						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
						if (comp == 3) {
						glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
						}
						else if (comp == 4) {
						glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
						}
						glBindTexture(GL_TEXTURE_2D, 0);
						stbi_image_free(image);*/
						textures.insert(std::make_pair(mp->diffuse_texname, bitmap));
					}
				}
			}
		}

		bmin[0] = bmin[1] = bmin[2] = std::numeric_limits<float>::max();
		bmax[0] = bmax[1] = bmax[2] = -std::numeric_limits<float>::max();

		{
			for (size_t s = 0; s < shapes.size(); s++) {
				DrawObject o;
				std::vector<float> vb;  // pos(3float), normal(3float), color(3float)
				for (size_t f = 0; f < shapes[s].mesh.indices.size() / 3; f++) {
					tinyobj::index_t idx0 = shapes[s].mesh.indices[3 * f + 0];
					tinyobj::index_t idx1 = shapes[s].mesh.indices[3 * f + 1];
					tinyobj::index_t idx2 = shapes[s].mesh.indices[3 * f + 2];

					int current_material_id = shapes[s].mesh.material_ids[f];

					if ((current_material_id < 0) || (current_material_id >= static_cast<int>(materials.size()))) {
						// Invaid material ID. Use default material.
						current_material_id = materials.size() - 1; // Default material is added to the last item in `materials`.
					}
					//if (current_material_id >= materials.size()) {
					//    std::cerr << "Invalid material index: " << current_material_id << std::endl;
					//}
					//
					float diffuse[3];
					for (size_t i = 0; i < 3; i++) {
						diffuse[i] = materials[current_material_id].diffuse[i];
					}
					float tc[3][2];
					if (attrib.texcoords.size() > 0) {
						assert(attrib.texcoords.size() > 2 * idx0.texcoord_index + 1);
						assert(attrib.texcoords.size() > 2 * idx1.texcoord_index + 1);
						assert(attrib.texcoords.size() > 2 * idx2.texcoord_index + 1);
						tc[0][0] = attrib.texcoords[2 * idx0.texcoord_index];
						tc[0][1] = 1.0f - attrib.texcoords[2 * idx0.texcoord_index + 1];
						tc[1][0] = attrib.texcoords[2 * idx1.texcoord_index];
						tc[1][1] = 1.0f - attrib.texcoords[2 * idx1.texcoord_index + 1];
						tc[2][0] = attrib.texcoords[2 * idx2.texcoord_index];
						tc[2][1] = 1.0f - attrib.texcoords[2 * idx2.texcoord_index + 1];
					}
					else {
						tc[0][0] = 0.0f;
						tc[0][1] = 0.0f;
						tc[1][0] = 0.0f;
						tc[1][1] = 0.0f;
						tc[2][0] = 0.0f;
						tc[2][1] = 0.0f;
					}

					float v[3][3];
					for (int k = 0; k < 3; k++) {
						int f0 = idx0.vertex_index;
						int f1 = idx1.vertex_index;
						int f2 = idx2.vertex_index;
						assert(f0 >= 0);
						assert(f1 >= 0);
						assert(f2 >= 0);
												
						v[0][k] = attrib.vertices[3 * f0 + k];
						v[1][k] = attrib.vertices[3 * f1 + k];
						v[2][k] = attrib.vertices[3 * f2 + k];
						//o.vertices.push_back(Vec3f(v[k][0], v[k][1], v[k][2]));
						bmin[k] = std::fmin(v[0][k], bmin[k]);
						bmin[k] = std::fmin(v[1][k], bmin[k]);
						bmin[k] = std::fmin(v[2][k], bmin[k]);
						bmax[k] = std::fmax(v[0][k], bmax[k]);
						bmax[k] = std::fmax(v[1][k], bmax[k]);
						bmax[k] = std::fmax(v[2][k], bmax[k]);
					}					

					float n[3][3];
					if (attrib.normals.size() > 0) {
						int f0 = idx0.normal_index;
						int f1 = idx1.normal_index;
						int f2 = idx2.normal_index;
						assert(f0 >= 0);
						assert(f1 >= 0);
						assert(f2 >= 0);
						for (int k = 0; k < 3; k++) {
							n[0][k] = attrib.normals[3 * f0 + k];
							n[1][k] = attrib.normals[3 * f1 + k];
							n[2][k] = attrib.normals[3 * f2 + k];
						}
					}
					else {
						// compute geometric normal
						CalcNormal(n[0], v[0], v[1], v[2]);
						n[1][0] = n[0][0];
						n[1][1] = n[0][1];
						n[1][2] = n[0][2];
						n[2][0] = n[0][0];
						n[2][1] = n[0][1];
						n[2][2] = n[0][2];
					}

					for (int k = 0; k < 3; k++) {
						vb.push_back(v[k][0]);
						vb.push_back(v[k][1]);
						vb.push_back(v[k][2]);						
						o.vertices.push_back(Vec3f(v[k][0], v[k][1], v[k][2]));
						vb.push_back(n[k][0]);
						vb.push_back(n[k][1]);
						vb.push_back(n[k][2]);
						// Combine normal and diffuse to get color.
						float normal_factor = 0.2;
						float diffuse_factor = 1 - normal_factor;
						float c[3] = {
							n[k][0] * normal_factor + diffuse[0] * diffuse_factor,
							n[k][1] * normal_factor + diffuse[1] * diffuse_factor,
							n[k][2] * normal_factor + diffuse[2] * diffuse_factor
						};
						float len2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
						if (len2 > 0.0f) {
							float len = sqrtf(len2);

							c[0] /= len;
							c[1] /= len;
							c[2] /= len;
						}
						vb.push_back(c[0] * 0.5 + 0.5);
						vb.push_back(c[1] * 0.5 + 0.5);
						vb.push_back(c[2] * 0.5 + 0.5);

						vb.push_back(tc[k][0]);
						vb.push_back(tc[k][1]);

						o.texcoords.push_back(Vec2f(tc[k][0], tc[k][1]));
					}
					
					o.indices.push_back(idx0.vertex_index);
					o.indices.push_back(idx1.vertex_index);
					o.indices.push_back(idx2.vertex_index);
				}

				//o.vb = 0;
				o.numTriangles = 0;

				// OpenGL viewer does not support texturing with per-face material.
				if (shapes[s].mesh.material_ids.size() > 0 && shapes[s].mesh.material_ids.size() > s) {
					// Base case
					o.material_id = shapes[s].mesh.material_ids[s];
				}
				else {
					o.material_id = materials.size() - 1; // = ID for default material.
				}

				if (vb.size() > 0) {
					/*glGenBuffers(1, &o.vb);
					glBindBuffer(GL_ARRAY_BUFFER, o.vb);
					glBufferData(GL_ARRAY_BUFFER, vb.size() * sizeof(float), &vb.at(0),
					GL_STATIC_DRAW);*/
					o.numTriangles = vb.size() / (3 + 3 + 3 + 2) / 3; // 3:vtx, 3:normal, 3:col, 2:texcoord

					printf("shape[%d] # of triangles = %d\n", static_cast<int>(s),
						o.numTriangles);
				}

				drawObjects->push_back(o);
			}
		}

		printf("bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
		printf("bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

		return true;
	}
	~ObjLoader();
};

#endif