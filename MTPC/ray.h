#pragma once
#include <glm/glm.hpp>
#include "geometry.h"
enum type {
	DIFFUSE, SPECULAR, TRANSMISSION
};

class Ray {
public:
	glm::vec3 start;
	glm::vec3 direction;
	enum type ray_type;
	Ray(glm::vec3 s, glm::vec3 d, type t) :start(s), direction(d), ray_type(t) {}
	Ray() {}
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
