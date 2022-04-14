duffman: duffman.c duffman.h
	cc -g -Wall $< -o $@

stress: stress.c duffman.h
	cc -g -Wall $< -o $@

test: stress
	./stress

dist: clean
	rm -f duffman.tar.gz && cd .. && tar czvf duffman/duffman.tar.gz duffman/*

clean:
	rm -f duffman *.o stress *.tar.gz *.asc
