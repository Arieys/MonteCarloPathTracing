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
//����
struct Vertex {
	// λ��
	glm::vec3 Position;
	// ������
	glm::vec3 Normal;
	// ��������
	glm::vec2 TexCoords;
	// u����
	glm::vec3 Tangent;
	// v����
	glm::vec3 Bitangent;
};
//����
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

class Face {//���Ӧ�Ķ�������
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
	/*  Model���� */
	//�洢��ĿǰΪֹ���ص����������Ż���ȷ�������ᱻ���ض�Ρ�
	vector<Texture> textures;
	std::vector <Face> faces;
	std::vector <Face> light;
	string directory;
	bool gammaCorrection;	
	Camera camera;

	/*  ����  */
	// ���캺��,��Ҫһ��3Dģ�͵��ļ�·��
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
	/*  ����   */
	// ���ļ�����֧��ASSIMP��չ��ģ�ͣ��������ɵ�����洢������ʸ���С�
	void loadModel(string const &path)
	{
		// ͨ��ASSIMP���ļ�
		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
		// ������
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) // �������0
		{
			cout << "����::ASSIMP:: " << importer.GetErrorString() << endl;
			return;
		}
		// �����ļ�·����Ŀ¼·��
		directory = path.substr(0, path.find_last_of('/'));

		// �Եݹ鷽ʽ����ASSIMP�ĸ��ڵ�
		processNode(scene->mRootNode, scene);
	}

	// �Եݹ鷽ʽ����ڵ㡣 ����λ�ڽڵ㴦��ÿ���������񣬲������ӽڵ㣨����У����ظ��˹��̡�
	void processNode(aiNode *node, const aiScene *scene)
	{
		// ����λ�ڵ�ǰ�ڵ��ÿ������
		for (unsigned int i = 0; i < node->mNumMeshes; i++)
		{
			// �ڵ��������������������������е�ʵ�ʶ���
			// ���������������ݣ��ڵ�ֻ��Ϊ������֯�ı��涫������ڵ�֮��Ĺ�ϵ����
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			processMesh(mesh, scene);
		}
		// �����Ǵ�����������������еĻ��������ǻ�ݹ鴦��ÿ���ӽڵ�
		for (unsigned int i = 0; i < node->mNumChildren; i++)
		{
			processNode(node->mChildren[i], scene);
		}
	}

	void processMesh(aiMesh *mesh, const aiScene *scene)
	{
		// ��ʱ����
		vector<Vertex> vertices;

		// ����ÿ������Ķ���
		for (unsigned int i = 0; i < mesh->mNumVertices; i++)
		{
			Vertex vertex;
			// ��������һ��ռλ����������Ϊassimpʹ�����Լ��������࣬����ֱ��ת��Ϊglm��vec3�࣬�����������Ƚ����ݴ��ݸ����ռλ��glm :: vec3��
			glm::vec3 vector;
			// λ��
			vector.x = mesh->mVertices[i].x;
			vector.y = mesh->mVertices[i].y;
			vector.z = mesh->mVertices[i].z;
			vertex.Position = vector;
			// ����
			vector.x = mesh->mNormals[i].x;
			vector.y = mesh->mNormals[i].y;
			vector.z = mesh->mNormals[i].z;
			vertex.Normal = vector;
			// ��������
			if (mesh->mTextureCoords[0]) // �����Ƿ�����������ꣿ
			{
				glm::vec2 vec;
				// �������ɰ���8����ͬ���������ꡣ ��ˣ����Ǽ������ǲ���ʹ�ö�����Ծ��ж�����������ģ�ͣ�����������ǲ��õ�һ�����ϣ�0����
				vec.x = mesh->mTextureCoords[0][i].x;
				vec.y = mesh->mTextureCoords[0][i].y;
				vertex.TexCoords = vec;
			}
			else
				vertex.TexCoords = glm::vec2(0.0f, 0.0f);
			// u����
			vector.x = mesh->mTangents[i].x;
			vector.y = mesh->mTangents[i].y;
			vector.z = mesh->mTangents[i].z;
			vertex.Tangent = vector;
			// v����
			vector.x = mesh->mBitangents[i].x;
			vector.y = mesh->mBitangents[i].y;
			vector.z = mesh->mBitangents[i].z;
			vertex.Bitangent = vector;
			vertices.push_back(vertex);
		}

		// �ӹ�����
		aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
		Material mat;
		aiColor3D color;
		double data_in;

		//��ȡmtl�ļ���������
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

		// 1. ��������ͼ
		int id = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");

		mat.Kd_id = id;

		//���ڱ���ÿ�������棨һ������һ�������ε����񣩲�������Ӧ�Ķ���������
		for (unsigned int i = 0; i < mesh->mNumFaces; i++)
		{
			aiFace face = mesh->mFaces[i];
			Face my_face;
			my_face.material = mat;
			// ����������������������Ǵ洢������������
			for (unsigned int j = 0; j < face.mNumIndices; j++)
				my_face.v.push_back(vertices[face.mIndices[j]]);
			glm::vec3 center;
			center = (my_face.v[0].Position + my_face.v[1].Position + my_face.v[2].Position) / 3.0f;
			my_face.morton_code = getMortonCode(center.x, center.y, center.z);
			my_face.calNormal();
			this->faces.push_back(my_face);
		}

		// ���Ǽ�����ɫ���еĲ���������Լ���� ÿ������������Ӧ����Ϊ'texture_diffuseN'������N�Ǵ�1��MAX_SAMPLER_NUMBER�����кš�
		//ͬ�����������������������ܽ᣺
		// diffuse: texture_diffuseN
		// specular: texture_specularN
		// normal: texture_normalN

		//// 2. �߹���ͼ
		//loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
		//// 3.������ͼ
		//loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
		//// 4. �߶���ͼ
		//loadMaterialTextures(material, aiTextureType_AMBIENT, "texture_height");

	}

	// ���������͵����в������������δ�����������������
	// ������Ϣ��ΪTexture�ṹ���ء�
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
			// ���֮ǰ�Ƿ��������������ǣ��������һ�ε�������������������
			bool skip = false;
			for (unsigned int j = 0; j < textures.size(); j++)
			{
				if (std::strcmp(textures[j].path.data(), str.C_Str()) == 0)
				{
					textures.push_back(textures[j]);
					skip = true;
					break;// �Ѽ��ؾ�����ͬ�ļ�·��������������һ�����Ż�����
				}
			}
			if (!skip)
			{   // �����δ���������������
				Texture texture;
				texture.img = cv::imread(directory + '/' + str.C_Str());
				texture.type = typeName;
				texture.path = str.C_Str();
				textures.push_back(texture);  //����洢Ϊ����ģ�ͼ��ص�������ȷ�����ǲ�������ظ�����
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
