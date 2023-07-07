#pragma once
#include <iostream>
#include <vector>
//#include <GL/glew.h>
//#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <string>
#include <fstream>
#include <map>
#include "morton code.h"
#include <opencv2/opencv.hpp>
//#include <atlimage.h>

//assimp
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

double dmin(double p1, double p2, double p3);
double dmax(double p1, double p2, double p3);
//顶点
struct Vertex {
	// 位置
	glm::vec3 Position;
	// 法向量
	glm::vec3 Normal;
	// 纹理坐标
	glm::vec2 TexCoords;
	// u向量
	glm::vec3 Tangent;
	// v向量
	glm::vec3 Bitangent;
};
//纹理
struct Texture {
	cv::Mat img;
	string type;
	string path;
};

struct Material
{
public:
	int Kd_id;
	glm::vec4 Ke;
	glm::vec4 Ka;
	glm::vec4 Kd;
	glm::vec4 Ks;
	glm::vec4 Tf;
	double Tr;
	double Ns;
	double Ni;
	void print();
};

class Face {//面对应的顶点索引
public:
	vector <Vertex> v;
	glm::vec3 normal;
	Material material;
	unsigned int morton_code;
	void print();
	double calAera();
	void calNormal();
};

class Camera {
public:
	glm::vec3 eye, look_at, up;
	double fovy;
	int width, height;
	void print();
};

class scene_data
{
public:
	/*  Model数据 */
	//存储到目前为止加载的所有纹理，优化以确保纹理不会被加载多次。
	vector<Texture> textures;
	std::vector <Face> faces;
	std::vector <Face> light;
	string directory;
	bool gammaCorrection;	
	Camera camera;

	/*  函数  */
	// 构造汉化,需要一个3D模型的文件路径
	scene_data(string const &path, bool gamma = false) : gammaCorrection(gamma)
	{
		read_scene(path);
	}

	bool read_scene(std::string filename) {
		loadModel(filename + ".obj");
		read_xml(filename + ".camera");
		handleLight();
		std::cout << "Face count = " << faces.size() << std::endl;
		std::cout << "Light face count = " << light.size() << std::endl;
		std::cout << "Texture count = " << textures.size() << std::endl;
	}
	bool read_xml(std::string filename);

	void handleLight()
	{
		for (int i = 0; i < faces.size(); i++) {
			if (faces[i].material.Ke != glm::vec4(0.0f, 0.0f, 0.0f, 1.0f)) {
				light.push_back(faces[i]);
			}
		}
	}

private:
	/*  函数   */
	// 从文件加载支持ASSIMP扩展的模型，并将生成的网格存储在网格矢量中。
	void loadModel(string const &path)
	{
		// 通过ASSIMP读文件
		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
		// 检查错误
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) // 如果不是0
		{
			cout << "错误::ASSIMP:: " << importer.GetErrorString() << endl;
			return;
		}
		// 检索文件路径的目录路径
		directory = path.substr(0, path.find_last_of('/'));

		// 以递归方式处理ASSIMP的根节点
		processNode(scene->mRootNode, scene);
	}

	// 以递归方式处理节点。 处理位于节点处的每个单独网格，并在其子节点（如果有）上重复此过程。
	void processNode(aiNode *node, const aiScene *scene)
	{
		// 处理位于当前节点的每个网格
		for (unsigned int i = 0; i < node->mNumMeshes; i++)
		{
			// 节点对象仅包含索引用来索引场景中的实际对象。
			// 场景包含所有数据，节点只是为了有组织的保存东西（如节点之间的关系）。
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			processMesh(mesh, scene);
		}
		// 在我们处理完所有网格（如果有的话）后，我们会递归处理每个子节点
		for (unsigned int i = 0; i < node->mNumChildren; i++)
		{
			processNode(node->mChildren[i], scene);
		}
	}

	void processMesh(aiMesh *mesh, const aiScene *scene)
	{
		// 临时数据
		vector<Vertex> vertices;

		// 遍历每个网格的顶点
		for (unsigned int i = 0; i < mesh->mNumVertices; i++)
		{
			Vertex vertex;
			// 我们声明一个占位符向量，因为assimp使用它自己的向量类，它不直接转换为glm的vec3类，所以我们首先将数据传递给这个占位符glm :: vec3。
			glm::vec3 vector;
			// 位置
			vector.x = mesh->mVertices[i].x;
			vector.y = mesh->mVertices[i].y;
			vector.z = mesh->mVertices[i].z;
			vertex.Position = vector;
			// 法线
			vector.x = mesh->mNormals[i].x;
			vector.y = mesh->mNormals[i].y;
			vector.z = mesh->mNormals[i].z;
			vertex.Normal = vector;
			// 纹理坐标
			if (mesh->mTextureCoords[0]) // 网格是否包含纹理坐标？
			{
				glm::vec2 vec;
				// 顶点最多可包含8个不同的纹理坐标。 因此，我们假设我们不会使用顶点可以具有多个纹理坐标的模型，因此我们总是采用第一个集合（0）。
				vec.x = mesh->mTextureCoords[0][i].x;
				vec.y = mesh->mTextureCoords[0][i].y;
				vertex.TexCoords = vec;
			}
			else
				vertex.TexCoords = glm::vec2(0.0f, 0.0f);
			// u向量
			vector.x = mesh->mTangents[i].x;
			vector.y = mesh->mTangents[i].y;
			vector.z = mesh->mTangents[i].z;
			vertex.Tangent = vector;
			// v向量
			vector.x = mesh->mBitangents[i].x;
			vector.y = mesh->mBitangents[i].y;
			vector.z = mesh->mBitangents[i].z;
			vertex.Bitangent = vector;
			vertices.push_back(vertex);
		}

		// 加工材料
		aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
		Material mat;
		aiColor3D color;
		double data_in;

		//读取mtl文件顶点数据
		material->Get(AI_MATKEY_COLOR_AMBIENT, color);
		mat.Ka = glm::vec4(color.r, color.g, color.b, 1.0);
		material->Get(AI_MATKEY_COLOR_DIFFUSE, color);
		mat.Kd = glm::vec4(color.r, color.g, color.b, 1.0);		
		material->Get(AI_MATKEY_COLOR_SPECULAR, color);
		mat.Ks = glm::vec4(color.r, color.g, color.b, 1.0);
		material->Get(AI_MATKEY_COLOR_TRANSPARENT, color);
		mat.Tf = glm::vec4(color.r, color.g, color.b, 1.0);
		material->Get(AI_MATKEY_COLOR_EMISSIVE, color);
		mat.Ke = glm::vec4(color.r, color.g, color.b, 1.0);
		cout << mat.Ke.x << " " << mat.Ke.y << " " << mat.Ke.z << std::endl;
		material->Get(AI_MATKEY_SHININESS, data_in);
		mat.Ns = data_in;
		material->Get(AI_MATKEY_REFRACTI, data_in);
		mat.Ni = data_in;
		material->Get(AI_MATKEY_OPACITY, data_in);
		mat.Tr = data_in;

		// 1. 漫反射贴图
		int id = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");

		mat.Kd_id = id;

		//现在遍历每个网格面（一个面是一个三角形的网格）并检索相应的顶点索引。
		for (unsigned int i = 0; i < mesh->mNumFaces; i++)
		{
			aiFace face = mesh->mFaces[i];
			Face my_face;
			my_face.material = mat;
			// 检索面的所有索引并将它们存储在索引向量中
			for (unsigned int j = 0; j < face.mNumIndices; j++)
				my_face.v.push_back(vertices[face.mIndices[j]]);
			glm::vec3 center;
			center = (my_face.v[0].Position + my_face.v[1].Position + my_face.v[2].Position) / 3.0f;
			my_face.morton_code = getMortonCode(center.x, center.y, center.z);
			my_face.calNormal();
			this->faces.push_back(my_face);
		}

		// 我们假设着色器中的采样器名称约定。 每个漫反射纹理应命名为'texture_diffuseN'，其中N是从1到MAX_SAMPLER_NUMBER的序列号。
		//同样适用于其他纹理，如下列总结：
		// diffuse: texture_diffuseN
		// specular: texture_specularN
		// normal: texture_normalN

		//// 2. 高光贴图
		//loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
		//// 3.法线贴图
		//loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
		//// 4. 高度贴图
		//loadMaterialTextures(material, aiTextureType_AMBIENT, "texture_height");

	}

	// 检查给定类型的所有材质纹理，如果尚未加载纹理，则加载纹理。
	// 所需信息作为Texture结构返回。
	int loadMaterialTextures(aiMaterial *mat, aiTextureType type, string typeName)
	{
		int material_cnt = mat->GetTextureCount(type);
		if (material_cnt > 1) 
			std::cout << "ERROR : MATERIAL CNT > 1" << std::endl;
		else if (material_cnt == 0) 
			return -1;
		for (unsigned int i = 0; i < mat->GetTextureCount(type); i++)
		{
			aiString str;
			mat->GetTexture(type, i, &str);
			// 检查之前是否加载了纹理，如果是，则继续下一次迭代：跳过加载新纹理
			bool skip = false;
			for (unsigned int j = 0; j < textures.size(); j++)
			{
				if (std::strcmp(textures[j].path.data(), str.C_Str()) == 0)
				{
					textures.push_back(textures[j]);
					skip = true;
					break;// 已加载具有相同文件路径的纹理，继续下一个（优化）。
				}
			}
			if (!skip)
			{   // 如果尚未加载纹理，请加载它
				Texture texture;
				texture.img = cv::imread(directory + '/' + str.C_Str());
				texture.type = typeName;
				texture.path = str.C_Str();
				textures.push_back(texture);  //将其存储为整个模型加载的纹理，以确保我们不会加载重复纹理。
			}
			return textures.size()-1;
		}

	}
};


class boundingBox {
public:
	double max_x, max_y, max_z, min_x, min_y, min_z; //boundingBox
	void print() {
		cout << min_x << " " << min_y << " " << min_z << " " << max_x << " " << max_y << " " << max_z << endl;
	}
};

bool compare(Face t1, Face t2);

enum type {
	DIFFUSE, SPECULAR, TRANSMISSION
};

class Ray {
public:
	glm::vec3 start;
	glm::vec3 direction;
	enum type ray_type;
	Ray(glm::vec3 s, glm::vec3 d, type t) :start(s), direction(d), ray_type(t){}
	Ray(){}
};

class image {
public:
	double *img;
	double *avg_img;
	int width, height;
	image(int width, int height) {
		this->width = width, this->height = height;
		img = new double[width * height * 3];
		avg_img = new double[width * height * 3];
		memset(img, 0.0, sizeof(double)*width*height * 3);//initialize img
		memset(avg_img, 0.0, sizeof(double)*width*height * 3);//initialize img
	}
	~image() {
		delete[]img;
		delete[]avg_img;
	}
	int getIndex(int x, int y, int rgb) { //x colunmn y row
		return y * width * 3 + 3 * x + rgb;
	}
	void average(int core) {
		for (int i = core; i + core < width; i++) {
			for (int j = core; j + core < height; j++) {
				double current_r = 0, current_g = 0, current_b = 0;
				for (int k = -core; k <= core; k++) {
					for (int m = -core; m <= core; m++) {
						current_r += img[getIndex(i + k, j + m, 0)];
						current_g += img[getIndex(i + k, j + m, 1)];
						current_b += img[getIndex(i + k, j + m, 2)];
					}
				}
				current_r /= pow(2 * core + 1, 2), current_g /= pow(2 * core + 1, 2), current_b /= pow(2 * core + 1, 2);
				avg_img[getIndex(i, j, 0)] = current_r;
				avg_img[getIndex(i, j, 1)] = current_g;
				avg_img[getIndex(i, j, 2)] = current_b;
			}
		}
	}
};

//record the intersection information
class intersection {
public:
	Ray ray;  //ray = start_point + t * direction
	glm::vec3 p; //intersection point
	glm::vec3 pn;//intersection point normal
	double t; //ray parameters t
	Face f;   //intersection face
};

//judge if a ray intersects with a triangle, ret is the intersection point
bool intersect(Ray &ray, Face &triangle, glm::vec3 &ret);

//judge if a ray intersects with a boungingbox
bool intersect(Ray &ray, boundingBox &b);
