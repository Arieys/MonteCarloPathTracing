#include <iostream>
#include "morton code.h"
#include <algorithm>
#include "svpng.inc"
#include <filesystem>
#include <chrono>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "pathTracer.h"

#include <glm/glm.hpp>

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

bool render_scene(std::string path, std::string filename, int N_ray_per_pixel)
{
	clock_t start,end;
	double duration;
	start = clock();
	scene_data scene(path + filename);
	sort(scene.faces.begin(), scene.faces.end(), compare);
	cout << "sort by morton code successfully" << endl;
	BVH bvh(scene);
	end = clock();
	std::cout << end << " " << start << std::endl;
	duration = static_cast<double>(end - start)/ CLOCKS_PER_SEC * 1000;

	std::cout << "Phase 1(read scene + bvh build) time cost = " << duration << " ms" << std::endl;
	
	image img(scene.camera.width, scene.camera.height); //定义图像空间

	start = clock();

	PathTracer pt;
	pt.generateImg(scene, bvh, img, N_ray_per_pixel);

	end = clock();

	duration = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000;

	std::cout << "Phase 2(ray tracing) = " << duration << " ms" << std::endl;

	pt.outputTime();
	std::string output_path = "../result/";
	imshow(img.img, scene.camera.width, scene.camera.height, output_path + filename, N_ray_per_pixel);

	return true;
}


int main()
{
	std::string path = std::filesystem::current_path().string() + "/../scene/";
	std::string filename1 = "bedroom";
	std::string filename2 = "veach-mis";
	std::string filename3 = "cornell-box";
	std::string filename4 = "scene01";
	std::cout << std::filesystem::current_path() << std::endl;
	render_scene(path, filename4, 100);
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

