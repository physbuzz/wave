CPP=g++

wave: wave.cpp ImageUtil.h utils.h
	$(CPP) -g wave.cpp -o wave

run: wave-fast
	./wave-fast

run-dbg: wave
	./wave

wave-fast: wave.cpp ImageUtil.h utils.h
	$(CPP) -O3 wave.cpp -o wave-fast

clean: 
	-rm wave

#tidy:
#	clang-tidy circlepack.cpp --
#
#clean: 
#	-rm ./circlepack
#	-rm ./fast
#
#FORCE: 
