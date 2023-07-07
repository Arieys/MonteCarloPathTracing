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
	}
	camera.print();
	return true;
}

void Material::print()
{
	std::cout << "Kd_id = " << Kd_id << std::endl;
	std::cout << "Kd: " << Kd.x << " " << Kd.y << " " << Kd.z << std::endl;
	std::cout << "Ks: " << Ks.x << " " << Ks.y << " " << Ks.z << std::endl;
	std::cout << "Ns: " << Ns << std::endl;
	std::cout << "Ni: " << Ni << std::endl;
}

void Face::print()
{
	for (int i = 0; i < this->v.size(); i++) {
		std::cout << "Position: " << v[i].Position.x << v[i].Position.y << v[i].Position.z << std::endl;
		std::cout << "Normal: " << v[i].Normal.x << v[i].Normal.y << v[i].Normal.z << std::endl;
	}
	std::cout << "morton code: " << morton_code << std::endl;
}

void Camera::print()
{
	std::cout << "eye: " << eye.x << " " << eye.y << " " << eye.z << std::endl;
	std::cout << "lookat: " << look_at.x << " " << look_at.y << " " << look_at.z << std::endl;
	std::cout << "up: " << up.x << " " << up.y << " " << up.z << std::endl;
	std::cout << "fovy: " << fovy << " width: " << width << " height: " << height << std::endl;
}

bool compare(Face t1, Face t2)
{
	return t1.morton_code < t2.morton_code;
}

bool intersect(Ray &ray, Face &triangle, glm::vec3 &ret)
{
	double t;
	//if (triangle.material[0] == 'P') std::cout << "ray direction: ", ray.direction.print();
	glm::vec3 norm = triangle.normal;
	//if (triangle.material[0] == 'P') std::cout <<"norm: ", norm.print();
	t = glm::dot(triangle.v[0].Position - ray.start, norm) / glm::dot(norm,ray.direction);

	glm::vec3 p = ray.start + static_cast<float>(t) * ray.direction;
	//if (triangle.material[0] == 'P') cout << "p: ", p.print();
	glm::vec3 ap = p - triangle.v[0].Position, bp = p - triangle.v[1].Position, cp = p - triangle.v[2].Position;
	glm::vec3 ab = triangle.v[1].Position - triangle.v[0].Position, bc = triangle.v[2].Position - triangle.v[1].Position, ca = triangle.v[0].Position - triangle.v[2].Position;
	glm::vec3 cross1 = glm::cross(ab,ap), cross2 = glm::cross(bc,bp), cross3 = glm::cross(ca,cp);
	double dir1 = glm::dot(cross1, norm), dir2 = glm::dot(cross2 ,norm), dir3 = glm::dot(cross3, norm);
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


double Face::calAera()
{
	double a = sqrt(glm::dot(v[1].Position - v[0].Position, v[1].Position - v[0].Position));
	double b = sqrt(glm::dot(v[2].Position - v[0].Position, v[2].Position - v[0].Position));
	double c = sqrt(glm::dot(v[2].Position - v[1].Position, v[2].Position - v[1].Position)); //length of a,b,c
	double cos_c = (a * a + b * b - c * c) / (2 * a * b); //ÓàÏÒ¶¨Àí
	double sin_c = sqrt(1 - pow(cos_c, 2));
	double aera = a * b * sin_c / 2;
	return aera;
}

void Face::calNormal()
{
	normal = glm::normalize(glm::cross(v[0].Position - v[1].Position, v[2].Position - v[1].Position));

}
