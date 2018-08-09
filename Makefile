CC= g++  -Wall
INCLUDE= -I /home/yerui/src/SeqLib/ -I /home/yerui/src/SeqLib/htslib/ /home/yerui/src/SeqLib/bin/libseqlib.a /home/yerui/src/SeqLib/bin/libhts.a /home/yerui/src/SeqLib/bin/libbwa.a /home/yerui/src/SeqLib/bin/libfml.a
LIBS= -llzma -lbz2 -L. -lgzstream -lz -lpthread

bamDCS:  main.cc main.h compass_search.o
	$(CC) $< compass_search.o $(INCLUDE) $(LIBS) -o $@
compass_search.o:

.PHONY: clean

clean:
	-rm bamDCS compass_search.o
