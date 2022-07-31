.PHONY: build

build: select-data merge-data cuberenderer


select-data: src/select-data.cpp src/csv_parser.hpp src/star_record.hpp
	g++ -o2 -g -o select-data src/select-data.cpp

merge-data: src/merge-data.cpp src/star_record.hpp
	g++ -o2 -g -o merge-data src/merge-data.cpp

cuberenderer: src/cuberenderer.cpp src/star_record.hpp src/ciedata.h
	g++ -o3 -g -mcmodel=medium -o cuberenderer src/cuberenderer.cpp -lIlmImf-2_5
