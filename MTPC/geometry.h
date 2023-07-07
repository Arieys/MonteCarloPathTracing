#pragma once
#include <string>
#include <vector>

#include "material.h"

#include <glm/glm.hpp>

#include <opencv2/opencv.hpp>
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
	std::string type;
	std::string path;
};

class Face {//���Ӧ�Ķ�������
public:
	std::vector <Vertex> v;
	glm::vec3 normal;
	Material material;
	unsigned int morton_code;
	double calAera();
	void calNormal();
};

bool compare(Face t1, Face t2);
double dmin(double p1, double p2, double p3);
double dmax(double p1, double p2, double p3);
