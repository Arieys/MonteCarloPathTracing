# MonteCarloPathTracing
该项目运行在VS2017环境下，使用前需要配置glm、eigen、svpng库，使用的依赖库版本为eigen-3.4.0，glm-0.9.9.8，glfw-3.3.6，glew-2.1.0  
在scene文件夹中放入需要渲染的场景文件，需要包含.obj，.mtl，.camera文件以及mtl文件中对应的纹理贴图 
  
其中.camera文件的格式如下  ：
eye 3.456 1.212 3.299 //相机位置  
lookat 2.699 1.195 2.645 //相机看向位置  
up -0.013 1.000 -0.011 //相机up方向  
fovy 39.4305 //相机fovy角度  
width 1280 //场景分辨率宽度  
height 720 //场景分辨率高度  
mtlname Light 16.4648 16.4648 16.4648 //光源radiance  
  
渲染生成的文件在.result文件夹中  
render_scene函数为渲染场景主函数，接收两个参数：文件名和每个像素的采样数，在调用时修改采样数以得到不同质量的光追结果

算法详解和代码详细解释见博客：https://blog.csdn.net/Listoree/article/details/124081645

PathTracer渲染结果展示
![bedroom-SPP256](https://user-images.githubusercontent.com/49404329/162612869-a46830ad-87bd-47f7-9ad7-1a823dc40e6b.png)  
bedroom-SPP256
![veach-mis-SPP100](https://user-images.githubusercontent.com/49404329/162612896-e3361e6f-d737-4975-aca9-82fcda2da1f8.png)  
veach-mis-SPP100
