# MonteCarloPathTracing
该项目运行在VS2017环境下，使用前需要配置glm、eigen、svpng库  
使用的依赖库版本为eigen-3.4.0，glm-0.9.9.8，glfw-3.3.6，glew-2.1.0  
配置完依赖库后在scene文件夹中放入需要渲染的场景文件，需要包含.obj，.mtl，.camera文件以及mtl文件中对应的纹理贴图，并在main函数里修改filename为场景文件名    
其中.camera文件的格式如下  
eye 3.456 1.212 3.299 //相机位置  
lookat 2.699 1.195 2.645 //相机看向位置  
up -0.013 1.000 -0.011 //相机up方向  
fovy 39.4305 //相机fovy角度  
width 1280 //场景分辨率宽度  
height 720 //场景分辨率高度  
mtlname Light 16.4648 16.4648 16.4648 //光源radiance    
渲染生成的文件在.result文件夹中  
