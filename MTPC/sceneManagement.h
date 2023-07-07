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

private:
	// 从文件加载支持ASSIMP扩展的模型，并将生成的网格存储在网格矢量中。
	void loadModel(string const& path);

	// 以递归方式处理节点。 处理位于节点处的每个单独网格，并在其子节点（如果有）上重复此过程。
	void processNode(aiNode* node, const aiScene* scene);

	void processMesh(aiMesh* mesh, const aiScene* scene);

	// 检查给定类型的所有材质纹理，如果尚未加载纹理，则加载纹理。
	int loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

	//读场景数据和相机数据
	bool read_scene(std::string filename);
	bool read_xml(std::string filename);

	//处理光源信息
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

