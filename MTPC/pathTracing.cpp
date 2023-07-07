#include "pathTracing.h"

bool russian_Roulette(double pr)
{
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double>u1(0, 1);
	double rnd = u1(e);
	//std::cout << rnd << std::endl;
	if (rnd < pr) return true;
	else return false;
}

bool Refract(glm::vec3 dir, glm::vec3 &normal, double eta, glm::vec3 &refract_dir) {
	// Ref: https://www.cnblogs.com/night-ride-depart/p/7429618.html
	glm::vec3 &i = dir;
	glm::vec3 &n = normal;
	float cosi = glm::dot(i, n);
	float cost2 = 1.0f - eta * eta * (1.0f - cosi * cosi);
	if (cost2 >= 0.0f) {
		refract_dir = i * static_cast<float>(eta) - n * static_cast<float>((eta * cosi + std::sqrt(cost2)));
		return true;
	}
	else {
		// total internal reflection
		return false;
	}
}

//importance sampling based on BRDF
glm::vec3 BRDFImportanceSampling(glm::vec3 &direction, type type, Face &f, double Ns) {
	// Ref : https://inst.eecs.berkeley.edu/~cs283/sp13/lectures/283-lecture11.pdf
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double>u1(0, 1);

	double phi = u1(e) * 2 * pi;
	double theta;
	if (type == DIFFUSE) { //diffuse
		//std::cout << "sample diffuse :" << f.material << std::endl;
		// for cosine-weighted Lambertian
		theta = asin(sqrt(u1(e)));
	}
	else if (type == SPECULAR) {
		// for sampling specular term
		//std::cout << "sample speculat :" << f.material << std::endl;
		theta = acos(pow(u1(e), (double)1 / (Ns + 1)));
	}
	else {
		cerr << "unknown sample type" << endl;
	}
	glm::vec3 sample(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi));
	glm::vec3 front;
	if (fabs(direction.x) > fabs(direction.y)) {
		front = glm::normalize(glm::vec3(direction.z, 0, -direction.x));
	}
	else {
		front = glm::normalize(glm::vec3(0, -direction.z, direction.y));
	}
	glm::vec3 right = glm::cross(direction,front);
	glm::vec3 ret = glm::normalize((right * sample.x) + (direction * sample.y) + (front * sample.z));
	return ret;
}

Ray nextRay(intersection &p, glm::vec3 dir,glm::vec3 kd, scene_data &data)
{
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double>u2(0, 1);

	Material *m = &p.f.material;

	glm::vec3 direction;

	//for refraction

	if (m->Ni > 1) {
		//std::cout << "fraction :" << m->name << std::endl;
		double n1, n2;
		double cos_in = glm::dot(-dir,p.pn);
		glm::vec3 normal;
		if (cos_in > 0) {
			// out of glass
			//std::cout << "out" << std::endl;
			normal = -p.pn;
			n1 = m->Ni;
			n2 = 1.0;
		}
		else {
			//std::cout << "in" << std::endl;
			// in to the glass
			normal = p.pn;
			n1 = 1.0;
			n2 = m->Ni;
		}

		// Ref : https://en.wikipedia.org/wiki/Schlick%27s_approximation
		double rf0 = pow((n1 - n2) / (n1 + n2),2);
		double fresnel = rf0 + (1.0f - rf0) * pow(1.0f - std::abs(cos_in),5);
		if (fresnel < u2(e) ){
			if (Refract(-dir, normal, n1 / n2, direction)) {
				return Ray(p.p, direction, TRANSMISSION);
			}
			else {
				// only specular for refraction material
				glm::vec3 reflect; //reflect ray direction
				glm::vec3 incoming = -dir;
				reflect = incoming - normal * glm::dot(incoming, normal) * 2.0f;
				return Ray(p.p, reflect, SPECULAR);
				// only specular for refraction material
			}

		}
	}

	double kd_norm = kd.length(), ks_norm = m->Ks.length();
	type ray_type;
	if (ks_norm != 0 && kd_norm / ks_norm < u2(e)) {
		//sample specular
		glm::vec3 reflect; //reflect ray direction
		glm::vec3 incoming = -dir;
		reflect = incoming - p.pn * glm::dot(incoming, p.pn) * 2.0f;
		direction = BRDFImportanceSampling(reflect, SPECULAR,p.f, m->Ns);
		ray_type = SPECULAR;
	}
	else {
		//sample diffuse
		direction = BRDFImportanceSampling(p.pn, DIFFUSE, p.f, m->Ns);
		ray_type = DIFFUSE;
	}

	Ray r(p.p + (direction * 0.01f), direction, ray_type);
	return r;
}

//calculate p -> dir radiance
glm::vec3 shade(intersection &p, glm::vec3 dir, scene_data &data, BVH &bvh)
{	
	glm::vec3 L_dir(0,0,0);		
	//处理自发光项（光源）
	if (p.f.material.Ke != glm::vec4(0,0,0,1.0)) {
		std::cout << "Hit light" << std::endl;
		L_dir = p.f.material.Ke;
		return L_dir;
	}

	//calculate the material infomation
	Material m = p.f.material;
	glm::vec3 kd;
	if (m.Kd_id != -1) {
		glm::vec3 garcov = findGarCor(p.f, p.p); //重心坐标
		//garcov.print();
		double row = p.f.v[0].TexCoords.x * garcov.x + p.f.v[1].TexCoords.x * garcov.y + p.f.v[2].TexCoords.x * garcov.z; //插值得到纹理坐标
		//std::cout << p.f.vt1.vtx << " " << p.f.vt2.vtx << " " << p.f.vt3.vtx << std::endl;
		double col = p.f.v[0].TexCoords.y * garcov.x + p.f.v[1].TexCoords.y * garcov.y + p.f.v[2].TexCoords.y * garcov.z;
		double irow = row - floor(row), icol = col - floor(col); //保留小数部分
		int r = irow * data.textures[m.Kd_id].img.size[1], c = icol * data.textures[m.Kd_id].img.size[0];
		//std::cout << m->name << " " << r << " " << c << std::endl;
		cv::Vec3b vec_3 = data.textures[m.Kd_id].img.at<cv::Vec3b>(r, c); //读取纹理值
		kd.x = (double)vec_3[2]/255, kd.y = (double)vec_3[1]/255, kd.z = (double)vec_3[0]/255; //得到Kd
	}
	else {
		kd = m.Kd;
	}

	//direct illumination
	//Uniformly sample light at x'
	
	L_dir.x = L_dir.y = L_dir.z = 0;
	static std::default_random_engine e(time(NULL));
	Face sample_face;
	for (int i = 0; i < data.light.size(); i++) {
		//for every direct light
		glm::vec3 xl; //xl is the sample point on light
		glm::vec3 vn; //normal at point xl
		Face l = data.light[i]; //l record the light attribute

		double total_aera = l.calAera();

		static std::uniform_real_distribution<double>u1(0, total_aera);
		double rnd = u1(e);
		//std::cout << rnd << std::endl;
		//uniformly generate a random number based on aera to choose which triagnle mesh to sample


		//sample a vertex on triangle mesh j
		sample_face = l;
		static std::uniform_real_distribution<double>u2(0, 1);
		double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
		float p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
		xl = sample_face.v[0].Position * p1 + sample_face.v[1].Position * p2 + sample_face.v[2].Position * p3;
		vn = sample_face.v[0].Normal * p1 + sample_face.v[1].Normal * p2 + sample_face.v[2].Normal * p3;

		glm::vec3 direction = glm::normalize(xl - p.p); //direction from point p to light

		//judge visibility
		float visibility = 1.0f;
		Ray rl;
		rl.start = p.p + (direction * 0.01f);//add micro turbulance
		rl.direction = direction;
		intersection inter;

		ray_intersect(rl, data, bvh, inter);
		if (inter.f.morton_code != sample_face.morton_code) {
			visibility = 0.0f;
		}

		if (glm::dot(direction, p.pn) > 0){
			float pdf_light = double(1)/total_aera; //pdf_light = 1/A where A is the aera of light source
			float cos_theta = abs(glm::dot(direction, vn) / direction.length() / vn.length());
			float cos_theta_hat = abs(glm::dot(direction, p.pn) / direction.length() / p.pn.length());
			float dist = max(1.0f,distance(xl, p.p));
			glm::vec3 intensity = l.material.Ke  *  cos_theta * cos_theta_hat / pow(dist, 2) / pdf_light * visibility;
			double kd_dots = glm::dot(direction,p.pn); //cos between light to intersection and face normal
			
			//only add diffuse
			if (kd_dots > 0) {
				L_dir.x += kd.x * intensity.x * kd_dots / pi;
				L_dir.y += kd.y * intensity.y * kd_dots / pi;
				L_dir.z += kd.z * intensity.z * kd_dots / pi;
			}
		}
	}

	//indirect illumination
	glm::vec3 L_indir(0, 0, 0);
	float P_RR = 0.6; //russian roulette probability

	//BRDF sample the hemisphere, get direction
	if (russian_Roulette(P_RR)) {
		Ray r = nextRay(p, dir, kd, data); //sample next ray based on BRDF
		//find the first object ray hits
		intersection ret;
		if (ray_intersect(r, data, bvh, ret)) {
			glm::vec3 intensity = shade(ret, -r.direction,data,bvh) / P_RR;

			if (r.ray_type == DIFFUSE) {
				if (ret.f.material.Ke == glm::vec4(0,0,0,1.0)) { //hit none emitting object
					L_indir.x += kd.x * intensity.x;
					L_indir.y += kd.y * intensity.y;
					L_indir.z += kd.z * intensity.z;
				}
			}
			else if (r.ray_type == SPECULAR) {
				L_indir.x += m.Ks.x *  intensity.x;
				L_indir.y += m.Ks.y *  intensity.y;
				L_indir.z += m.Ks.z *  intensity.z;
			}
			else {
				L_indir.x += m.Tf.x * intensity.x;
				L_indir.y += m.Tf.y * intensity.y;
				L_indir.z += m.Tf.z * intensity.z;
			}
		}
	}

	return L_dir + L_indir;
}

//clock_t start_intersect, end_intersect;
//clock_t start_shade, end_shade;
//double intersect_duration = 0, shade_duration = 0;
clock_t s, e;
double duration = 0;
//逐像素发射光线
void generateImg(scene_data &scene, BVH &bvh, image &img, int N_ray_per_pixel)
{
	scene.camera.up = glm::normalize(scene.camera.up);
	glm::vec3 dir = scene.camera.look_at - scene.camera.eye;
	float l = dir.length();
	float dy = tan(scene.camera.fovy/2/180 * pi) * l; //根据fovy求世界坐标系中半屏幕的x,y大小
	float dx = dy / scene.camera.height * scene.camera.width;

	std::cout << l << " " << dy << " " << dx << std::endl;

	glm::vec3 screen_center = scene.camera.look_at; //屏幕在世界坐标系下中心点的坐标
	float pdx = 2 * dx / scene.camera.width, pdy = 2 * dy / scene.camera.height; //每个像素在世界坐标系下的大小
	
	glm::vec3 screen_x_dir = glm::normalize(glm::cross(dir, scene.camera.up));

	glm::vec3 screen_y_dir = scene.camera.up;

	glm::vec3 screen_pdy = screen_y_dir * pdy;
	glm::vec3 screen_pdx = screen_x_dir * pdx;

	glm::vec3 start_point = screen_center - (screen_x_dir * dx) + (scene.camera.up * dy); //逐像素发射光线 像素的中心位置
	glm::vec3 pos = start_point;
	for (int i = 0; i < scene.camera.height; i++) {
		pos = start_point - (screen_pdy * static_cast<float>(i));
		
		for (int j = 0; j < scene.camera.width; j++) {
			mutex m;
			glm::vec3 current_radiance(0, 0, 0);
			//mutex clock_inte_m, clock_shade_m;
			omp_set_num_threads(min(N_ray_per_pixel,8));
			#pragma omp parallel for
			for (int k = 0; k < N_ray_per_pixel; k++) {
				Ray ray;
				ray.start = scene.camera.eye, ray.direction = glm::normalize(pos - scene.camera.eye);
				//generate ray
			
				intersection ret;
				if (ray_intersect(ray, scene, bvh, ret)) {
					glm::vec3 radiance = shade(ret, -ray.direction, scene, bvh);

					m.lock();
					current_radiance.x += radiance.x/N_ray_per_pixel;
					current_radiance.y += radiance.y/N_ray_per_pixel;
					current_radiance.z += radiance.z/N_ray_per_pixel;
					m.unlock();
				}
			}
			img.img[img.getIndex(j, i, 0)] = current_radiance.x;
			img.img[img.getIndex(j, i, 1)] = current_radiance.y;
			img.img[img.getIndex(j, i, 2)] = current_radiance.z;

			pos = pos + screen_pdx;
		}
		std::cout << "img pos: (" << i << ",x) generate successfully " << std::endl;
	}
	std::cout << "Garvo time cost = " << duration << " ms" << std::endl;
}

//ray intersection detect using bvh
void bvh_intersect(Ray ray, BVH &bvh, intersection &v, int current_node, int current_level, bool &flag)
{
	int index = bvh.findIndex(current_node, current_level);
	if (intersect(ray, bvh.bvh[index].b)) {
		if (bvh.bvh[index].level == bvh.Level) {
				//reach leaf node
				glm::vec3 ret;
				if (intersect(ray, *bvh.bvh[index].object, ret)) {
					//if(bvh.bvh[index].object->material[0] == 'P')
						//std::cout << "hit " << bvh.bvh[index].object->material << std::endl;
					intersection i;
					i.p = ret;
					i.ray = ray;
					i.t = (i.p.x - ray.start.x) / ray.direction.x;
					i.f = *bvh.bvh[index].object;

					glm::vec3 garcov = findGarCor(i.f, i.p);
					i.pn = (i.f.v[0].Normal * garcov.x) + (i.f.v[1].Normal * garcov.y) + (i.f.v[2].Normal * garcov.z);
					//if (ray.direction * i.f.norm < 0) i.pn = i.f.norm;
						//else i.pn = i.f.norm.negative();

					if (flag == false) {
						if (i.t > 0) {
							v = i;
							flag = true;
						}
					}
					else if (i.t > 0 && i.t < v.t) v = i; //remain the smallest t one
				}
				//}

				return;
			}	
		else {
				int next_node1 = 2 * current_node + 1;
				int next_node2 = 2 * current_node + 2;
				bvh_intersect(ray, bvh, v, next_node1, current_level+1,flag);
				bvh_intersect(ray, bvh, v, next_node2, current_level+1,flag);
		}
	}
}

bool comp(intersection &i1, intersection &i2)
{
	return i1.t < i2.t;
}

//ray intersect with scene
bool ray_intersect(Ray ray, scene_data &scene, BVH &bvh, intersection &ret)
{
	//start_intersect = clock();
	bool flag = false;
	bvh_intersect(ray, bvh, ret, 0, 0, flag);
	//end_intersect = clock();
	//intersect_duration += static_cast<double>(end_intersect - start_intersect) / CLOCKS_PER_SEC * 1000;
	return flag;
}


//求重心坐标
glm::vec3 findGarCor(Face &f, glm::vec3 p)
{
	s = clock();
	//Plan 1 works but time consuming, use eigen to sovle a linear function
	//Eigen::Matrix <double, 4, 3> A;
	//Eigen::Matrix <double, 4, 1> B;
	//Eigen::MatrixXd res;


	//A << f.v1.x, f.v2.x, f.v3.x,
	//	f.v1.y, f.v2.y, f.v3.y,
	//	f.v1.z, f.v2.z, f.v3.z,
	//	1, 1, 1;
	//B << p.x, p.y, p.z, 1;

	//res = A.colPivHouseholderQr().solve(B);

	//Vertex ret;
	//ret.x = res(0, 0), ret.y = res(1, 0), ret.z = res(2, 0);


	//plan 2 faster than plan 1 but no idea how it works
	glm::vec3 e1, e2, e3, d1, d2, d3;
	e1 = f.v[2].Position - f.v[1].Position;
	e2 = f.v[0].Position - f.v[2].Position;
	e3 = f.v[1].Position - f.v[0].Position;
	d1 = p - f.v[0].Position, d2 = p - f.v[1].Position, d3 = p - f.v[2].Position;
	glm::vec3 n = glm::cross(e1,e2);
	double an = glm::dot(n,n);
	double b1, b2, b3;
	b1 = glm::dot(glm::cross(e1, d3), n) / an;
	b2 = glm::dot(glm::cross(e2, d1), n) / an;
	b3 = glm::dot(glm::cross(e3, d2), n) / an;
	glm::vec3 ret;
	ret.x = b1, ret.y = b2, ret.z = b3;
	e = clock();
	duration += static_cast<double>(e - s) / CLOCKS_PER_SEC * 1000;
	return ret;
}