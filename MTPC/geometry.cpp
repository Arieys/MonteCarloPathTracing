#include "geometry.h"

double Face::calAera()
{
	double a = sqrt(glm::dot(v[1].Position - v[0].Position, v[1].Position - v[0].Position));
	double b = sqrt(glm::dot(v[2].Position - v[0].Position, v[2].Position - v[0].Position));
	double c = sqrt(glm::dot(v[2].Position - v[1].Position, v[2].Position - v[1].Position)); //length of a,b,c
	double cos_c = (a * a + b * b - c * c) / (2 * a * b); //”‡œ“∂®¿Ì
	double sin_c = sqrt(1 - pow(cos_c, 2));
	double aera = a * b * sin_c / 2;
	return aera;
}

void Face::calNormal()
{
	normal = glm::normalize(glm::cross(v[0].Position - v[1].Position, v[2].Position - v[1].Position));

}

bool compare(Face t1, Face t2)
{
	return t1.morton_code < t2.morton_code;
}

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