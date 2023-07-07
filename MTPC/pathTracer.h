#pragma once
#include "sceneManagement.h"
#include "BVH.h"
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include <omp.h>

#define MULTI_THREAD

const double pi = 3.1415926;

class PathTracer {
public:
	PathTracer() { triangle_intersection_time = bvh_intersection_time = shading_time = garcor_time = 0; }
	void generateImg(scene_data &scene, BVH &bvh, image &img, int N_ray_per_pixel);
	void setMaxDepth(int m) { this->MAXDEPTH = m; }
	void outputTime();
private:
	mutex timer_mutex;
	double triangle_intersection_time;
	double bvh_intersection_time;
	double shading_time;
	double garcor_time;

	int MAXDEPTH = 6;

	//test russian_roulette with probability pr
	bool russian_Roulette(double pr);

	//calculate p -> w0 radiance
	glm::vec3 shade(intersection &p, glm::vec3 dir, scene_data &data, BVH &bvh, int depth);

	//shading related
	bool Refract(glm::vec3 dir, glm::vec3& normal, double eta, glm::vec3& refract_dir);
	glm::vec3 BRDFImportanceSampling(glm::vec3& direction, type type, Face& f, double Ns);
	Ray nextRay(intersection& p, glm::vec3 dir, glm::vec3 kd, scene_data& data);

	//geometry trlated
	glm::vec3 findGarCor(Face &f, glm::vec3 p);	
	bool comp(intersection &i1, intersection &i2);
	//judge if a ray intersects with a triangle, ret is the intersection point
	bool intersect(Ray& ray, Face& triangle, glm::vec3& ret);
	//judge if a ray intersects with a boungingbox
	bool intersect(Ray& ray, boundingBox& b);	
	void bvh_intersect(Ray ray, BVH &bvh, intersection &v, int current_node, int current_level, bool &flag);
	bool ray_intersect(Ray ray, scene_data &scene, BVH &bvh, intersection &ret);
};
