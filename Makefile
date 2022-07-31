.PHONY: build

build: select-data merge-data cuberenderer


select-data: src/select-data.cpp src/csv_parser.hpp src/star_record.hpp
	g++ -o2 -g -o select-data src/select-data.cpp

merge-data: src/merge-data.cpp src/star_record.hpp
	g++ -o2 -g -o merge-data src/merge-data.cpp

cuberenderer: src/cuberenderer.cpp src/star_record.hpp src/ciedata.h debug/recordsize src/temperature_lookup.hpp
	g++ -O3 -march=native -fprofile-dir=. -fprofile-generate -fallow-store-data-races -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fcx-limited-range -fexcess-precision=fast -g -mcmodel=medium -o debug/cuberenderer-profile src/cuberenderer.cpp -lIlmImf-2_5
	[ ! -f 'debug/#root#stars#cuberenderer.gcda' ] || rm 'debug/#root#stars#cuberenderer.gcda'
	(cd debug && lzop -d <../stars.lzop | dd bs="`./recordsize`0000" count=16384 | ./cuberenderer-profile)
	g++ -O3 -march=native -fprofile-dir=debug -fprofile-use -fallow-store-data-races -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fcx-limited-range -fexcess-precision=fast -g -mcmodel=medium -o cuberenderer src/cuberenderer.cpp -lIlmImf-2_5

debug/recordsize: src/star_record.hpp
	mkdir -p debug
	g++ -g -o debug/recordsize src/recordsize.cpp

debug/gen-temperature-lookup: src/gen-temperature-lookup.cpp src/star_record.hpp
	g++ -O3 -march=native -fallow-store-data-races -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fcx-limited-range -fexcess-precision=fast -g -o debug/gen-temperature-lookup src/gen-temperature-lookup.cpp

src/temperature_lookup.hpp: debug/gen-temperature-lookup stars.lzop
	pv stars.lzop | lzop -d | debug/gen-temperature-lookup > src/temperature_lookup.hpp
