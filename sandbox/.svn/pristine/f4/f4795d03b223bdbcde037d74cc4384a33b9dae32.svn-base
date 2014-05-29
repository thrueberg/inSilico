solution "Beta Cells"
	configurations { "Debug", "Local", "Cluster" }

project "Beta Cells"
	language		 "C++" 
	kind			 "ConsoleApp" 
	files			{ "**.hpp","**.cpp" }
	links		 	{ "gomp" }
	defines			{"NTHREADS=1"}
	buildoptions 	{
					  "-std=c++11",
					  "-fopenmp",
					  "-I $$INSILICOROOT"
					}


	configuration 	 	{ "Debug" }
		defines 	 	{ "_DEBUG", "DEBUG" }
		flags		 	{ "Symbols" } 
		buildoptions	{ 	
							"-isystem /usr/include/eigen3" 
						}
		targetname	  	 "d_Beta_Cells" 


	configuration		{ "Local" }
		defines			{ "NDEBUG" }
		flags			{ "Optimize" }
		buildoptions	{ 	
							"-isystem /usr/include/eigen3" 
						}
		targetname		 "Beta_Cells"


	configuration		{ "Cluster" } 
		defines			{ "NDEBUG" }
		flags			{ "OptimizeSpeed" }
		buildoptions	{ 
							"-isystem ~/lib/Eigen"
						}
		targetname		 "Beta_Cells"

