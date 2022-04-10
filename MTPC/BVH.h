#pragma once
#include <vector>
#include "sceneManagement.h"
#include <algorithm>
using namespace std;


class BVHNode {
public:
	Face *object;
	boundingBox b;
	BVHNode(){}	
	void findBondingBox(Face &data);
	void findBondingBox(BVHNode &n1, BVHNode &v2);
	void findBondingBox(BVHNode &n1);
	double findMax(double p1, double p2, double p3, double p4);
	double findMin(double p1, double p2, double p3, double p4);
	int level;
};


class BVH{
public:
	BVH(scene_data &data);
	BVHNode* bvh;
	void buildBVH(vector <Face> &data);
	bool boundingBoxOverlap(int node1, int node2, int l);
	bool haveRightSubtree(int bvh, int l);
	int t, Nv, Nr, Nc, Lc, Lr, Lv, Level;
	int findIndex(int i, int l);
	void print() {
		std::cout << Level << std::endl;
		for (int i = 0; i < Nr; i++) {
			bool flag = true;
			if (bvh[i].object == NULL) flag = false;
			//std::cout << "Level: " <<  bvh[i].level << " ObjN: " << flag << " ";
			//bvh[i].b.print();
			if (bvh[i].b.min_y < 0) {
				std::cout << "Level = " << bvh[i].level << "  ";
				bvh[i].b.print();
			}
		}
	}
	~BVH()
	{
		delete[] bvh;
	}
};