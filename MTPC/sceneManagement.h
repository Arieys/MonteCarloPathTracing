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

#include "geometry.h"
#include "material.h"
#include "ray.h"
#include "boundingbox.h"
#include "camera.h"
//assimp
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

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

private:
	// ���ļ�����֧��ASSIMP��չ��ģ�ͣ��������ɵ�����洢������ʸ���С�
	void loadModel(string const& path);

	// �Եݹ鷽ʽ����ڵ㡣 ����λ�ڽڵ㴦��ÿ���������񣬲������ӽڵ㣨����У����ظ��˹��̡�
	void processNode(aiNode* node, const aiScene* scene);

	void processMesh(aiMesh* mesh, const aiScene* scene);

	// ���������͵����в������������δ�����������������
	int loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

	//���������ݺ��������
	bool read_scene(std::string filename);
	bool read_xml(std::string filename);

	//�����Դ��Ϣ
	void handleLight();
};

//rendered image class
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
	void average(int core);
};

