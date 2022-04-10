#include "pathTracing.h"

bool russian_Roulette(double pr)
{
	static std::default_random_engine e;
	static std::uniform_real_distribution<double>u1(0, 1);
	double rnd = u1(e);
	//std::cout << rnd << std::endl;
	if (rnd < pr) return true;
	else return false;
}

Vertex BRDFImportanceSampling(Vertex &direction, int type, Face &f, int Ns) {
	// Ref : https://inst.eecs.berkeley.edu/~cs283/sp13/lectures/283-lecture11.pdf
	static std::default_random_engine e;
	static std::uniform_real_distribution<double>u1(0, 1);

	float phi = u1(e) * 2 * pi;
	float theta;
	if (type == DIFFUSE) { //diffuse
		// for cosine-weighted Lambertian
		theta = asin(sqrt(u1(e)));
	}
	else if (type == SPECULAR) {
		// for sampling specular term
		theta = acos(pow(u1(e), 1 / (Ns + 1)));
	}
	else {
		cerr << "unknown sample type" << endl;
	}
	Vertex sample(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi));
	Vertex front;
	if (fabs(direction.x) > fabs(direction.y)) {
		front = Vertex(direction.z, 0, -direction.x);
		front.normalize();
	}
	else {
		front = Vertex(0, -direction.z, direction.y);
		front.normalize();
	}
	Vertex right = direction.cross(front);
	Vertex ret = (right * sample.x) + (direction * sample.y) + (front * sample.z);
	ret.normalize();
	return ret;
}


//calculate p -> dir radiance
Vertex shade(intersection &p, Vertex dir, scene_data &data, BVH &bvh)
{	
	Vertex L_dir(0,0,0);		
	//处理自发光项（光源）
	if (data.light_map.find(p.f.material) != data.light_map.end()) {
		L_dir = data.light_map[p.f.material]->radiance;
		return L_dir;
	}

	//calculate the material infomation
	Material *m = data.material_map[p.f.material];
	Vertex kd;
	if (m->mapping_flag == true) {
		Vertex garcov = findGarCor(p.f, p.p); //重心坐标
		//garcov.print();
		double row = p.f.vt1.vtx * garcov.x + p.f.vt2.vtx * garcov.y + p.f.vt3.vtx * garcov.z; //插值得到纹理坐标
		//std::cout << p.f.vt1.vtx << " " << p.f.vt2.vtx << " " << p.f.vt3.vtx << std::endl;
		double col = p.f.vt1.vty * garcov.x + p.f.vt2.vty * garcov.y + p.f.vt3.vty * garcov.z;
		double irow = row - floor(row), icol = col - floor(col); //保留小数部分
		int r = irow * m->map_height, c = icol * m->map_width;
		//std::cout << m->name << " " << r << " " << c << std::endl;
		cv::Vec3b vec_3 = m->img.at<cv::Vec3b>(r, c); //读取纹理值
		kd.x = (double)vec_3[2]/255, kd.y = (double)vec_3[1]/255, kd.z = (double)vec_3[0]/255; //得到Kd
	}
	else {
		kd = m->kd;
	}

	//direct illumination
	//Uniformly sample light at x'
	
	L_dir.x = L_dir.y = L_dir.z = 0;
	static std::default_random_engine e;
	Face sample_face;
	for (int i = 0; i < data.l.size(); i++) {
		//for every direct light
		Vertex xl; //xl is the sample point on light
		Vertex vn; //normal at point xl
		Light l = data.l[i]; //l record the light attribute

		double total_aera = 0;
		int n_triangle_face = data.material_map[data.l[i].name]->f.size();
		double *triangle_aera = new double[n_triangle_face];
		for (int j = 0; j < n_triangle_face; j++) {
			//for every triangle mesh
			total_aera += data.material_map[data.l[i].name]->f[j].calAera(); //this is for pdf_light
			triangle_aera[j] = total_aera;
		}
		static std::uniform_real_distribution<double>u1(0, total_aera);
		double rnd = u1(e);
		//std::cout << rnd << std::endl;
		//uniformly generate a random number based on aera to choose which triagnle mesh to sample
		for (int j = 0; j < n_triangle_face; j++) {
			if (rnd < triangle_aera[j]) {
				//sample a vertex on triangle mesh j
				sample_face = data.material_map[data.l[i].name]->f[j];
				static std::uniform_real_distribution<double>u2(0, 1);
				double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
				double p1 = rnd1 / (rnd1 + rnd2 + rnd3), p2 = rnd2 / (rnd1 + rnd2 + rnd3), p3 = rnd3 / (rnd1 + rnd2 + rnd3);
				xl = sample_face.v1 * p1 + sample_face.v2 * p2 + sample_face.v3 * p3;
				vn = sample_face.vn1 * p1 + sample_face.vn2 * p2 + sample_face.vn3 * p3;
				break;
			}
		}
		delete[]triangle_aera;
		Vertex direction = xl - p.p; //direction from point p to light
		direction.normalize();

		//判断是否遮挡
		double visibility = 1;
		Ray rl;
		rl.start = p.p + (direction * 0.01);
		rl.direction = direction;
		//rl.start.print(), rl.direction.print();
		intersection inter;
		ray_intersect(rl, data, bvh, inter);
		if (inter.f.material != sample_face.material) {
			//std::cout << inter.f.material << " " << p.f.material << std::endl;
			//inter.p.print();
			visibility = 0;
		}

		if (direction * p.pn > 0){
			double pdf_light = double(1)/total_aera; //pdf_light = 1/A where A is the aera of light source
			double cos_theta = abs(direction * vn / direction.norm() / vn.norm());
			double cos_theta_hat = abs(direction * p.pn / direction.norm() / p.pn.norm());
			Vertex intensity = l.radiance  *  cos_theta * cos_theta_hat / pow(distance(xl,p.p),2) / pdf_light * visibility;
			double kd_dots = direction * p.pn; //cos between light to intersection and face normal
			
			if (kd_dots > 0) {
				L_dir.x += kd.x * intensity.x * kd_dots;
				L_dir.y += kd.y * intensity.y * kd_dots;
				L_dir.z += kd.z * intensity.z * kd_dots;
			}
			Vertex incoming = direction.negative(); //光源到该点的光线方向
			Vertex reflect = incoming - p.pn * (incoming * p.pn) * 2; //反射光线的方向
			double ks_dots = reflect * dir; //与法向量夹角的cos
			if (ks_dots > 0) {
				L_dir.x += m->ks.x * pow(ks_dots, m->Ns) * intensity.x;
				L_dir.y += m->ks.y * pow(ks_dots, m->Ns) * intensity.y;
				L_dir.z += m->ks.z * pow(ks_dots, m->Ns) * intensity.z;
			}
		}
	}

	//return L_dir;
	//indirect illumination
	Vertex L_indir(0, 0, 0);
	double P_RR = 0.6; //russian roulette probability

	//Uniformly sample the hemisphere, get direction
	if (russian_Roulette(P_RR)) {
		Vertex hemi_direction;
		static std::uniform_real_distribution<double>u2(0, 2 * pi);
		static std::uniform_real_distribution<double>u3(0, pi);
		while (1) {
			double theta = u2(e), fai = u3(e); 
			hemi_direction.x = sin(theta) * cos(fai);
			hemi_direction.y = sin(theta) * sin(fai);
			hemi_direction.z = cos(theta);
			if (hemi_direction * p.pn > 0) break; //采样接受
		}

		double cos_theta2 = hemi_direction * p.pn;
		//std::cout << cos_theta2 << std::endl;
		double pdf_hemi = double(1) / (2 * pi);
		Ray r;
		r.start = p.p + (hemi_direction * 0.01);
		r.direction = hemi_direction;
		r.direction.normalize();
		//find the first object ray hits
		intersection ret;
		if (ray_intersect(r, data, bvh, ret)) {
			//std::cout << ret.f.material << " " << p.f.material << std::endl;
			if (data.light_map.find(ret.f.material) == data.light_map.end()) { //hit none emitting object
				Vertex intensity = shade(ret, hemi_direction.negative(),data,bvh) * cos_theta2 / pdf_hemi / P_RR;
				//intensity.print();
				double kd_dots = hemi_direction * p.pn; //cos between light to intersection and face normal
				if (kd_dots > 0) {
					L_indir.x += kd.x * intensity.x * kd_dots;
					L_indir.y += kd.y * intensity.y * kd_dots;
					L_indir.z += kd.z * intensity.z * kd_dots;
				}
				Vertex reflect; //reflect ray direction
				Vertex incoming = hemi_direction.negative();
				reflect = incoming - p.pn * (incoming * p.pn) * 2;
				double ks_dots = reflect * dir;
				//std::cout << ks_dots << std::endl;
				if (ks_dots > 0) {
					L_indir.x += m->ks.x * pow(ks_dots, m->Ns) * intensity.x * (1 + m->Ns);
					L_indir.y += m->ks.y * pow(ks_dots, m->Ns) * intensity.y * (1 + m->Ns);
					L_indir.z += m->ks.z * pow(ks_dots, m->Ns) * intensity.z * (1 + m->Ns);
				}
			}
			else {
				//std::cout << "hit light" << std::endl;
			}
		}
	}
	//if(L_dir.x > 1) std::cout << "Dir radiance: ", L_dir.print();
	//if(L_indir.x > 1) std::cout << "inDir radiance: ", L_indir.print();
	return L_dir + L_indir;
}

//逐像素发射光线
void generateImg(scene_data &scene, BVH &bvh, image &img)
{
	scene.camera.up.normalize();
	Vertex dir = scene.camera.look_at - scene.camera.eye;
	double l = dir.norm();
	double dy = tan(scene.camera.fovy/2/180 * pi) * l; //根据fovy求世界坐标系中半屏幕的x,y大小
	double dx = dy / scene.camera.height * scene.camera.width;

	std::cout << l << " " << dy << " " << dx << std::endl;

	Vertex screen_center = scene.camera.look_at; //屏幕在世界坐标系下中心点的坐标
	double pdx = 2 * dx / scene.camera.width, pdy = 2 * dy / scene.camera.height; //每个像素在世界坐标系下的大小
	
	Vertex screen_x_dir = dir.cross(scene.camera.up);
	screen_x_dir.normalize();
	Vertex screen_y_dir = scene.camera.up;

	Vertex screen_pdy = screen_y_dir * pdy;
	Vertex screen_pdx = screen_x_dir * pdx;

	Vertex start_point = screen_center - (screen_x_dir * dx) + (scene.camera.up * dy); //逐像素发射光线 像素的中心位置
	Vertex pos = start_point;
	for (int i = 0; i < scene.camera.height; i++) {
		pos = start_point - (screen_pdy * i);
		
		for (int j = 0; j < scene.camera.width; j++) {
			omp_set_num_threads(N_ray_per_pixel);
			#pragma omp parallel for
			for (int k = 0; k < N_ray_per_pixel; k++) {
				Ray ray;
				ray.start = pos, ray.direction = pos - scene.camera.eye;
				ray.direction.normalize();
				//generate ray
			
				intersection ret;
				if (ray_intersect(ray, scene, bvh, ret)) {
					Vertex radiance = shade(ret, ray.direction.negative(), scene, bvh);
					//if(radiance.x != 0) radiance.print();
					
					//radiance.print();

					if (radiance.x > 1) radiance.x = 1;
					if (radiance.y > 1) radiance.y = 1;
					if (radiance.z > 1) radiance.z = 1;

					img.img[img.getIndex(j, i, 0)] += radiance.x/N_ray_per_pixel;
					img.img[img.getIndex(j, i, 1)] += radiance.y/N_ray_per_pixel;
					img.img[img.getIndex(j, i, 2)] += radiance.z/N_ray_per_pixel;
				}
			}
			//std::cout << "img pos: (" << i << "," << j << ") generate successfully " << std::endl;
			pos = pos + screen_pdx;
		}
		std::cout << "img pos: (" << i << ",x) generate successfully " << std::endl;
	}
}

//ray intersection detect using bvh
void bvh_intersect(Ray ray, BVH &bvh, intersection &v, int current_node, int current_level, bool &flag)
{
	int index = bvh.findIndex(current_node, current_level);
	if (bvh.bvh[index].level == bvh.Level) {
		//reach leaf node
		Vertex ret;
		if (intersect(ray, *bvh.bvh[index].object, ret)) {
			intersection i;
			i.p = ret;
			i.ray = ray;
			i.t = (i.p.x - ray.start.x) / ray.direction.x;
			i.f = *bvh.bvh[index].object;
			if (ray.direction * i.f.norm < 0) i.pn = i.f.norm;
			else i.pn = i.f.norm.negative();
			if (flag == false) {
				if (i.t > 0) {
					v = i;
					flag = true;
				}
			}
			else if (i.t > 0 && i.t < v.t) v = i; //remain the smallest t one
		}
		return;
	}
	else {
		if (intersect(ray, bvh.bvh[index].b)) {
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
	bool flag = false;
	bvh_intersect(ray, bvh, ret, 0, 0, flag);
	return flag;
}


//求重心坐标
Vertex findGarCor(Face &f, Vertex p)
{
	Eigen::Matrix <double, 4, 3> A;
	Eigen::Matrix <double, 4, 1> B;
	Eigen::MatrixXd res;


	A << f.v1.x, f.v2.x, f.v3.x,
		f.v1.y, f.v2.y, f.v3.y,
		f.v1.z, f.v2.z, f.v3.z,
		1, 1, 1;
	B << p.x, p.y, p.z, 1;

	res = A.colPivHouseholderQr().solve(B);

	Vertex ret;
	ret.x = res(0, 0), ret.y = res(1, 0), ret.z = res(2, 0);

	return ret;
}