#pragma once
#include "sceneManagement.h"
#include "BVH.h"
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include <omp.h>

const double pi = 3.1415926;

//test russian_roulette with probability pr
bool russian_Roulette(double pr);

//calculate p -> w0 radiance
glm::vec3 shade(intersection &p, glm::vec3 dir, scene_data &data, BVH &bvh);

void bvh_intersect(Ray ray, BVH &bvh, intersection &v, int current_node, int current_level, bool &flag);

bool ray_intersect(Ray ray, scene_data &scene, BVH &bvh, intersection &ret);

bool comp(intersection &i1, intersection &i2);

void generateImg(scene_data &scene, BVH &bvh, image &img, int N_ray_per_pixel);

bool ray_intersect(Ray ray, scene_data &scene, BVH &bvh, intersection &ret);

glm::vec3 findGarCor(Face &f, glm::vec3 p);
