#include <iostream>
#include "sceneManagement.h"
#include "morton code.h"
#include <algorithm>
#include "pathTracing.h"
#include "svpng.inc"

// 输出 SRC 数组中的数据到图像
void imshow(double* SRC, int WIDTH, int HEIGHT, string filename, int N_ray_per_pixel)
{

	unsigned char* image = new unsigned char[WIDTH * HEIGHT * 3];// 图像buffer
	unsigned char* p = image;
	double* S = SRC;    // 源数据

	FILE* fp;
	char *buffer = new char[10];
	_itoa(N_ray_per_pixel, buffer, 10);
	fopen_s(&fp, (filename + "-SPP" + buffer + ".png").c_str(), "wb");

	for (int i = 0; i < HEIGHT; i++)
	{
		for (int j = 0; j < WIDTH; j++)
		{
			*p++ = (unsigned char)glm::clamp((*S++) * 255, 0.0, 255.0);  // R 通道
			*p++ = (unsigned char)glm::clamp((*S++) * 255, 0.0, 255.0);  // G 通道
			*p++ = (unsigned char)glm::clamp((*S++) * 255, 0.0, 255.0);  // B 通道
		}
	}

	svpng(fp, WIDTH, HEIGHT, image, 0);
}

bool render_scene(std::string filename, int N_ray_per_pixel)
{
	scene_data scene;
	
	scene.read_scene(filename);

	sort(scene.f.begin(), scene.f.end(), compare);
	cout << "sort by morton code successfully" << endl;
	BVH bvh(scene);
	
	image img(scene.camera.width, scene.camera.height); //定义图像空间

	generateImg(scene, bvh, img, N_ray_per_pixel);

	imshow(img.img, scene.camera.width, scene.camera.height, filename, N_ray_per_pixel);
}


int main()
{
	std::string filename1 = "bedroom";
	std::string filename2 = "veach-mis";
	std::string filename3 = "cornell-box";

	render_scene(filename3, 100);
	//render_scene(filename2, 10);
	//render_scene(filename3, 10);
	//render_scene(filename1, 100);
	//render_scene(filename2, 100);
	//render_scene(filename3, 100);
	//render_scene(filename1, 256);
	//render_scene(filename2, 256);
	//render_scene(filename3, 256);
	//render_scene(filename1, 1000);
	//render_scene(filename2, 1000);
	//render_scene(filename3, 1000);

}

