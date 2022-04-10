#include "sceneManagement.h"

double dmin(double p1, double p2, double p3)
{
	if (p1 <= p2 && p1 <= p3) return p1;
	else if (p2 <= p1 && p2 <= p3) return p2;
	else return p3;
}

double dmax(double p1, double p2, double p3)
{
	if (p1 >= p2 && p1 >= p3) return p1;
	else if (p2 >= p1 && p2 >= p3) return p2;
	else return p3;
}

bool scene_data::read_mtl(std::string filename)
{
	std::ifstream fin(filename);
	std::string readline;
	Material *pm = NULL;
	while (getline(fin, readline)) {
		if (strcmp(readline.substr(0, 6).c_str(), "newmtl") == 0) {
			//std::cout << readline.substr(7).c_str() << std::endl;
			pm = new Material();
			pm->name = readline.substr(7);
			pm->mapping_flag = false;
			material_map[pm->name] = pm; //set map
		}
		else if (strcmp(readline.substr(0, 2).c_str(), "Kd") == 0) {
			readline = readline.substr(3);
			int blank;
			blank = readline.find(" ");
			pm->kd.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			pm->kd.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			pm->kd.z = atof(readline.c_str());
			//std::cout << "Kd: " << pm->kd1 << " " << pm->kd2 << " " << pm->kd3 << std::endl;
		}
		else if (strcmp(readline.substr(0, 2).c_str(), "Ks") == 0) {
			readline = readline.substr(3);
			int blank;
			blank = readline.find(" ");
			pm->ks.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			pm->ks.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			pm->ks.z = atof(readline.c_str());
			//std::cout << "Ks: " << pm->ks1 << " " << pm->ks2 << " " << pm->ks3 << std::endl;
		}
		else if (strcmp(readline.substr(0, 2).c_str(), "Ns") == 0) {
			readline = readline.substr(3);
			pm->Ns = atof(readline.c_str());
			//std::cout << "Ns: " << pm->Ns << std::endl;
		}
		else if (strcmp(readline.substr(0, 2).c_str(), "Ni") == 0) {
			readline = readline.substr(3);
			pm->Ni = atof(readline.c_str());
			//std::cout << "Ni: " << pm->Ni << std::endl;
		}
		else if (strcmp(readline.substr(0, 6).c_str(), "map_Kd") == 0) {
			readline = readline.substr(7);
			pm->map_Kd_filename = readline;
			pm->mapping_flag = true;
			pm->readinMap(); //read in texture map
			//std::cout << "map_Kd: " << pm->map_Kd_filename << std::endl;
		}
	}
	std::cout << "number of materials = " << material_map.size() << std::endl;
	return true;
}

bool scene_data::read_obj(std::string filename)
{
	std::ifstream fin(filename);
	std::string readline;
	std::string material;
	while (getline(fin, readline)) {
		if (readline[0] == 'v' && readline[1] == ' ') {
			//read vertex coordinate info
			Vertex current_v;
			double v1, v2, v3;
			int blank;
			readline = readline.substr(2);
			blank = readline.find(" ");
			v1 = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			v2 = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			v3 = atof(readline.c_str());
			//std::cout << "v: " << v1 << " " << v2 << " " << v3 << std::endl;
			current_v.x = v1, current_v.y = v2, current_v.z = v3;
			v.push_back(current_v);
		}
		else if (readline[0] == 'v' && readline[1] == 'n' && readline[2] == ' ') {
			//read vertex normal info
			Vertex current_vn;
			double vn1, vn2, vn3;
			int blank;
			readline = readline.substr(3);
			blank = readline.find(" ");
			vn1 = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			vn2 = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			vn3 = atof(readline.c_str());
			//std::cout << "vn: " << vn1 << " " << vn2 << " " << vn3 << std::endl;
			current_vn.x = vn1, current_vn.y = vn2, current_vn.z = vn3;
			vn.push_back(current_vn);
		}
		else if(readline[0] == 'v' && readline[1] == 't' && readline[2] == ' '){
			//read vertex t info
			VertexT current_vt;
			double vt1, vt2;
			int blank;
			readline = readline.substr(3);
			blank = readline.find(" ");
			vt1 = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			vt2 = atof(readline.substr(0, blank).c_str());
			//std::cout << "vt: " << vt1 << " " << vt2 << std::endl;
			current_vt.vtx = vt1, current_vt.vty = vt2;
			vt.push_back(current_vt);
		}
		else if (strcmp(readline.substr(0, 6).c_str(), "usemtl") == 0) {
			material = readline.substr(7);
		}
		else if (readline[0] == 'f' && readline[1] == ' ') {
			//read face info
			int index1, index2, index3;
			Face f;
			f.material = material;
			int slash, blank;
			readline = readline.substr(2);
			slash = readline.find('/');
			index1 = atoi(readline.substr(0, slash).c_str()) - 1;
			f.v1 = v[atoi(readline.substr(0, slash).c_str()) - 1];
			readline = readline.substr(slash + 1);
			slash = readline.find('/');
			f.vn1 = vn[atoi(readline.substr(0, slash).c_str())-1];
			readline = readline.substr(slash + 1);
			blank = readline.find(' ');
			f.vt1 = vt[atoi(readline.substr(0, blank).c_str())-1];
			
			readline = readline.substr(blank + 1);
			slash = readline.find('/');
			index2 = atoi(readline.substr(0, slash).c_str()) - 1;
			f.v2 = v[atoi(readline.substr(0, slash).c_str())-1];
			readline = readline.substr(slash + 1);
			slash = readline.find('/');
			f.vn2 = vn[atoi(readline.substr(0, slash).c_str())-1];
			readline = readline.substr(slash + 1);
			blank = readline.find(' ');
			f.vt2 = vt[atoi(readline.substr(0, blank).c_str())-1];

			readline = readline.substr(blank + 1);
			slash = readline.find('/');
			index3 = atoi(readline.substr(0, slash).c_str()) - 1;
			f.v3 = v[atoi(readline.substr(0, slash).c_str())-1];
			readline = readline.substr(slash + 1);
			slash = readline.find('/');
			f.vn3 = vn[atoi(readline.substr(0, slash).c_str())-1];
			readline = readline.substr(slash + 1);
			f.vt3 = vt[atoi(readline.substr(0, slash).c_str())-1];
			//std::cout << index1 << " " << index2 << " " << index3 << std::endl;
			f.calNorm();
			//if (vn[f.v1.vn_index].vn.x != vn[f.v2.vn_index].vn.x) {
			//	std::cout << "normal unmatched: " << f.v1.v_index << " " << f.v2.v_index << std::endl;
			//}
			Vertex center;
			center = (f.v1 + f.v2 + f.v3) / 3;
			//center.print();
			f.morton_code = getMortonCode(center.x, center.y, center.z);
			//std::cout << f.morton_code << std::endl;
			this->f.push_back(f);
			material_map[material]->f.push_back(f);
		}
	}

	std::cout << "number of vertices = " << v.size() << std::endl;
	std::cout << "number of faces = " << f.size() << std::endl;
	return true;
}

bool scene_data::read_xml(std::string filename)
{
	std::ifstream fin(filename);
	std::string readline;
	while (getline(fin, readline)) {
		if (strcmp(readline.substr(0, 3).c_str(), "eye") == 0) {
			readline = readline.substr(4);
			int blank;
			blank = readline.find(" ");
			camera.eye.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			camera.eye.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			camera.eye.z = atof(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 6).c_str(), "lookat") == 0) {
			readline = readline.substr(7);
			int blank;
			blank = readline.find(" ");
			camera.look_at.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			camera.look_at.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			camera.look_at.z = atof(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 2).c_str(), "up") == 0) {
			readline = readline.substr(3);
			int blank;
			blank = readline.find(" ");
			camera.up.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			camera.up.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			camera.up.z = atof(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 4).c_str(), "fovy") == 0) {
			readline = readline.substr(5);
			camera.fovy = atof(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 5).c_str(), "width") == 0) {
			readline = readline.substr(6);
			camera.width = atoi(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 6).c_str(), "height") == 0) {
			readline = readline.substr(7);
			camera.height = atoi(readline.c_str());
		}
		else if (strcmp(readline.substr(0, 7).c_str(), "mtlname") == 0) {
			readline = readline.substr(8);
			Light *l = new Light;
			int blank;
			blank = readline.find(" ");
			l->name = readline.substr(0, blank);
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			l->radiance.x = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			blank = readline.find(" ");
			l->radiance.y = atof(readline.substr(0, blank).c_str());
			readline = readline.substr(blank + 1);
			l->radiance.z = atof(readline.c_str());
			this->l.push_back(*l);
			light_map[l->name] = l;
			l->print();
		}
	}
	camera.print();
	return true;
}

bool scene_data::read_scene(std::string filename)
{
	std::string objfilename = filename + ".obj";
	std::string mtlfilename = filename + ".mtl";
	std::string xmlfilename = filename + ".camera";
	this->read_mtl(mtlfilename);
	this->read_obj(objfilename);
	this->read_xml(xmlfilename);
	return true;
}

void Material::print()
{
	Material *pm = this;
	std::cout << this->name << std::endl;
	std::cout << "Kd: ";
	this->kd.print();
	std::cout << "Ks: ";
	this->ks.print();
	std::cout << "Ns: " << pm->Ns << std::endl;
	std::cout << "Ni: " << pm->Ni << std::endl;
	std::cout << "map_Kd: " << pm->map_Kd_filename << std::endl;
}

void Face::print()
{
	std::cout << "material: " << this->material << std::endl;
	v1.print(), v2.print(), v3.print();
	vn1.print(), vn2.print(), vn3.print();
	vt1.print(), vt2.print(), vt3.print();
	std::cout << "morton code: " << morton_code << std::endl;
}

void Camera::print()
{
	std::cout << "eye: " << eye.x << " " << eye.y << " " << eye.z << std::endl;
	std::cout << "lookat: " << look_at.x << " " << look_at.y << " " << look_at.z << std::endl;
	std::cout << "up: " << up.x << " " << up.y << " " << up.z << std::endl;
	std::cout << "fovy: " << fovy << " width: " << width << " height: " << height << std::endl;
}

void Light::print()
{
	std::cout << name << " radiance: " << radiance.x << " " << radiance.y << " " << radiance.z << std::endl;
}

bool compare(Face t1, Face t2)
{
	return t1.morton_code < t2.morton_code;
}

bool intersect(Ray &ray, Face &triangle, Vertex &ret)
{
	double t;
	//if (triangle.material[0] == 'P') std::cout << "ray direction: ", ray.direction.print();
	Vertex norm = triangle.norm;
	//if (triangle.material[0] == 'P') std::cout <<"norm: ", norm.print();
	t = ((triangle.v1 - ray.start) * norm) / (norm * ray.direction);

	Vertex p = ray.start + (ray.direction * t);
	//if (triangle.material[0] == 'P') cout << "p: ", p.print();
	Vertex ap = p - triangle.v1, bp = p - triangle.v2, cp = p - triangle.v3;
	Vertex ab = triangle.v2 - triangle.v1, bc = triangle.v3 - triangle.v2, ca = triangle.v1 - triangle.v3;
	Vertex cross1 = ab.cross(ap), cross2 = bc.cross(bp), cross3 = ca.cross(cp);
	double dir1 = cross1 * norm, dir2 = cross2 * norm, dir3 = cross3 * norm;
	//if(triangle.material[0] == 'P') std::cout << dir1 << " " << dir2 << " " << dir3 << std::endl;
	double j1, j2, j3;
	j1 = dir1 * dir2, j2 = dir1 * dir3, j3 = dir2 * dir3;
	//std::cout << j1 << " " << j2 << " " << j2;
	bool judge = false;
	if (j1 >= 0 && j2 >= 0 && j3 >= 0) judge = true;
	ret = p;
	return judge;
}

bool intersect(Ray &ray, boundingBox &b)
{
	//if ray parellel with axis
	//if (ray.direction.x == 0 && ray.direction.y == 0 && ray.direction.z == 0)return false;
	//if (ray.direction.x == 0 && ray.direction.y == 0) {
	//	// ray // z axis
	//	if (ray.start.z > b.max_z && ray.direction.z > 0) return false;
	//	if (ray.start.z < b.min_z && ray.direction.z < 0) return false;
	//	if (ray.start.x < b.min_x || ray.start.x > b.max_x || ray.start.y < b.min_y || ray.start.y > b.max_y) return false;
	//	return true;
	//}
	//if (ray.direction.x == 0 && ray.direction.z == 0) {
	//	// ray // y axis
	//	if (ray.start.y > b.max_y && ray.direction.y > 0) return false;
	//	if (ray.start.y < b.min_y && ray.direction.y < 0) return false;
	//	if (ray.start.x < b.min_x || ray.start.x > b.max_x || ray.start.z < b.min_z || ray.start.z > b.max_z) return false;
	//	return true;
	//}
	//if (ray.direction.z == 0 && ray.direction.y == 0) {
	//	// ray // x axis
	//	if (ray.start.y > b.max_y && ray.direction.y > 0) return false;
	//	if (ray.start.y < b.min_y && ray.direction.y < 0) return false;
	//	if (ray.start.z < b.min_z || ray.start.z > b.max_z || ray.start.y < b.min_y || ray.start.y > b.max_y) return false;
	//	return true;
	//}
	//if ray not parellel with axis
	//ray.direction x y z != 0
	double txmin, txmax, tymin, tymax, tzmin, tzmax;
	txmin = (b.min_x - ray.start.x) / ray.direction.x;
	txmax = (b.max_x - ray.start.x) / ray.direction.x;
	tymin = (b.min_y - ray.start.y) / ray.direction.y;
	tymax = (b.max_y - ray.start.y) / ray.direction.y;
	tzmin = (b.min_z - ray.start.z) / ray.direction.z;
	tzmax = (b.max_z - ray.start.z) / ray.direction.z;
	if (txmin > txmax) {
		double temp = txmin;
		txmin = txmax, txmax = temp;
	}
	if (tymin > tymax) {
		double temp = tymin;
		tymin = tymax, tymax = temp;
	}
	if (tzmin > tzmax) {
		double temp = tzmin;
		tzmin = tzmax, tzmax = temp;
	}

	if (txmax < 0 || tymax < 0 || tzmax < 0) return false;
	if (txmin <= 0 && tymin <= 0 && tzmin <= 0) return true;
	if (dmax(txmin, tymin, tzmin) <= dmin(txmax, tymax, tzmax)) return true;
	else return false;
}

//calculate the distance between v1 and v2
double distance(Vertex v1, Vertex v2)
{
	return (v1 - v2).norm();
}

double Face::calAera()
{
	double a = (v2 - v1).norm(), b = (v3 - v1).norm(), c = (v3-v2).norm(); //length of a,b,c
	double cos_c = (a * a + b * b - c * c) / (2 * a * b); //ÓàÏÒ¶¨Àí
	double sin_c = sqrt(1 - pow(cos_c, 2));
	double aera = a * b * sin_c / 2;
	return aera;
}

void Face::calNorm()
{
	norm = (v1 - v2).cross(v3 - v1);
	norm.normalize();
}
