#include "BVH.h"

int count_level_20 = 0;

int count_set_bits(int x)
{
	int res = 0;
	for (int i = 0; i < 32; i++) {
		if (x & 0x1 == 1) res++;
		x >>= 1;
	}
	return res;
}


double BVHNode::findMax(double p1, double p2, double p3, double p4)
{

	if (p1 > p2 && p1 > p3 && p1 > p4) return p1;
	if (p2 > p1 && p2 > p3 && p2 > p4) return p2;
	if (p3 > p1 && p3 > p2 && p3 > p4) return p3;
	return p4;
}

double BVHNode::findMin(double p1, double p2, double p3, double p4)
{

	if (p1 < p2 && p1 < p3 && p1 < p4) return p1;
	if (p2 < p1 && p2 < p3 && p2 < p4) return p2;
	if (p3 < p1 && p3 < p2 && p3 < p4) return p3;
	return p4;
}


BVH::BVH(scene_data &data)
{
	buildBVH(data.faces);
	cout << "Build BVH success" << endl;
}

void BVH::buildBVH(vector <Face> &data)
{
	t = data.size();  //总object数
	Lc = pow(2, ceil(log2(data.size())));             //叶节点的个数
	Lv = Lc - t;                                      //虚叶节点的个数
	Nc = 2 * Lc - 1;//总节点数
	Nv = 2 * Lv - count_set_bits(Lv);//总虚节点数
	Nr = 2 * t - 1 + count_set_bits(Lv); //总实节点数
	Level = floor(log2(Nc));   //最大层数
	cout << "Total real = " << Nr << endl;
	bvh = new BVHNode[Nr];

	for (int l = Level; l >= 0; l--) {
		int current_level_v = Lv >> (Level - l);
		int start = pow(2, l) - 1, end = pow(2, l + 1) - 1 - current_level_v;
		int k = 0;
		cout << start << " " << end-1 << " " << endl;
		for (int i = start; i < end; i++) {
			if (l == Level) {
				//叶节点
				bvh[findIndex(i, l)].object = &data[k];
				bvh[findIndex(i, l)].findBondingBox(data[k]);
				bvh[findIndex(i, l)].level = l;
				k++;
			}
			else {
				bvh[findIndex(i, l)].object = NULL;
				if (haveRightSubtree(i,l)) {
					bvh[findIndex(i, l)].findBondingBox(bvh[findIndex(2 * i + 1, l + 1)], bvh[findIndex(2 * i + 2, l + 1)]);
					bvh[findIndex(i, l)].level = l;
				}
				else {
					//右子树缺失
					bvh[findIndex(i, l)].findBondingBox(bvh[findIndex(2 * i + 1, l + 1)]);
					bvh[findIndex(i, l)].level = l;
				}
				
			}
		}
		cout << "level " << l << "build successfully" << endl;
	}
}

void BVHNode::findBondingBox(Face &data)
{
	glm::vec3 p1 = data.v[0].Position, p2 = data.v[1].Position, p3 = data.v[2].Position;
	
	b.max_x = dmax(p1.x, p2.x, p3.x);
	b.max_y = dmax(p1.y, p2.y, p3.y);
	b.max_z = dmax(p1.z, p2.z, p3.z);
	b.min_x = dmin(p1.x, p2.x, p3.x);
	b.min_y = dmin(p1.y, p2.y, p3.y);
	b.min_z = dmin(p1.z, p2.z, p3.z);
}

int BVH::findIndex(int i, int l)
{
	int Lvl = Lv >> (Level - l + 1);
	int Nvl = 2 * Lvl - count_set_bits(Lvl);
	return i - Nvl;
}

void BVHNode::findBondingBox(BVHNode &v1, BVHNode &v2)
{
	this->b.max_x = max(v1.b.max_x, v2.b.max_x);
	this->b.max_y = max(v1.b.max_y, v2.b.max_y);
	this->b.max_z = max(v1.b.max_z, v2.b.max_z);
	this->b.min_x = min(v1.b.min_x, v2.b.min_x);
	this->b.min_y = min(v1.b.min_y, v2.b.min_y);
	this->b.min_z = min(v1.b.min_z, v2.b.min_z);
}

void BVHNode::findBondingBox(BVHNode &v1)
{
	this->b.max_x = v1.b.max_x;
	this->b.max_y = v1.b.max_y;
	this->b.max_z = v1.b.max_z;
	this->b.min_x = v1.b.min_x;
	this->b.min_y = v1.b.min_y;
	this->b.min_z = v1.b.min_z;
}

bool BVH::haveRightSubtree(int bvh, int l)
{
	if (2 * bvh + 2 >= pow(2, l + 2) - 1 - (Lv >> (Level-l- 1))) {
		return false;
	}
	else return true;
}

bool BVH::boundingBoxOverlap(int node1, int node2,int l)
{
	double minx, miny, minz, maxx, maxy, maxz;
	minx = max(bvh[findIndex(node1, l)].b.min_x, bvh[findIndex(node2, l)].b.min_x);
	miny = max(bvh[findIndex(node1, l)].b.min_y, bvh[findIndex(node2, l)].b.min_y);
	minz = max(bvh[findIndex(node1, l)].b.min_z, bvh[findIndex(node2, l)].b.min_z);
	maxx = min(bvh[findIndex(node1, l)].b.max_x, bvh[findIndex(node2, l)].b.max_x);
	maxy = min(bvh[findIndex(node1, l)].b.max_y, bvh[findIndex(node2, l)].b.max_y);
	maxz = min(bvh[findIndex(node1, l)].b.max_z, bvh[findIndex(node2, l)].b.max_z);
	if (minx > maxx || miny > maxy || minz > maxz) return false;
	else return true;
}