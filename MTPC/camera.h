#pragma once
#include <iostream>

#include <glm/glm.hpp>

class Camera {
public:
	glm::vec3 eye, look_at, up;
	double fovy;
	int width, height;
	void print();
};