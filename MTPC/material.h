#pragma once
#include <glm/glm.hpp>
struct Material
{
public:
	int Kd_id;
	glm::vec4 Ke;
	glm::vec4 Ka;
	glm::vec4 Kd;
	glm::vec4 Ks;
	double Tr;
	double Ns;
	double Ni;
};
