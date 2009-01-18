.phony: tags

all: Makefile
	make -f Makefile
Makefile: CMakeLists.txt
	cmake -DITK_DIR=$(ITK_DIR) -DCMAKE_BUILD_TYPE=Debug .
tags:
	ctags --langmap=c++:+.txx --languages=c++ -R $(ITK_SOURCE) CImg.h


