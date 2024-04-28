CPP=g++

wave: wave.cpp ImageUtil.h utils.h
	$(CPP) -g wave.cpp -o wave

run: wave-fast
	./wave-fast

run-dbg: wave
	./wave

wave-fast: wave.cpp ImageUtil.h utils.h
	$(CPP) -O3 wave.cpp -o wave-fast -fopenmp

clean: 
	-rm wave

# make wave-fast && mv wave-fast vid && cd vid && rm ./*.bmp && ./wave-fast && ffmpeg -r 60 -i out%04d.bmp -vcodec libx264 -pix_fmt yuv420p wave-newtonian.mp4
#
#tidy:
#	clang-tidy circlepack.cpp --
#
#clean: 
#	-rm ./circlepack
#	-rm ./fast
#
#FORCE: 
