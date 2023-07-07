#pragma once
#include <iostream>

class boundingBox {
public:
	double max_x, max_y, max_z, min_x, min_y, min_z; //boundingBox
	void print() {
		std::cout << min_x << " " << min_y << " " << min_z << " " << max_x << " " << max_y << " " << max_z << std::endl;
	}
};
