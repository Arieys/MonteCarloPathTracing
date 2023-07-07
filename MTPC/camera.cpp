#include "camera.h"

void Camera::print()
{
	std::cout << "eye: " << eye.x << " " << eye.y << " " << eye.z << std::endl;
	std::cout << "lookat: " << look_at.x << " " << look_at.y << " " << look_at.z << std::endl;
	std::cout << "up: " << up.x << " " << up.y << " " << up.z << std::endl;
	std::cout << "fovy: " << fovy << " width: " << width << " height: " << height << std::endl;
}