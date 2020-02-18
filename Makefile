CC= g++ -O3 -Wall
INCLUDE= -I /home/yerui/src/SeqLib/ -I /home/yerui/src/SeqLib/htslib/ /home/yerui/src/SeqLib/bin/libseqlib.a /home/yerui/src/SeqLib/bin/libhts.a /home/yerui/src/SeqLib/bin/libbwa.a /home/yerui/src/SeqLib/bin/libfml.a -I /home/yerui/anaconda3/include -I /home/yerui/src/alglib-3.14.0/src -L /home/yerui/anaconda3/lib -L. /home/yerui/src/alglib-3.14.0/src/*.o
LIBS= -llzma -lbz2 -L. -lgzstream -lz -lpthread

bamDCS:  main.cc main.h
	$(CC) $< -o $@ $(INCLUDE) $(LIBS)
compass_search.o:

.PHONY: clean

clean:
	-rm bamDCS
