all: libfoo.so libbar.so

libfoo.so: foo.c
	$(CC) $(CFLAGS) -shared -fPIC -o libfoo.so foo.c

libbar.so: bar.c
	$(CC) $(CFLAGS) -shared -fPIC -o libbar.so bar.c

tidy:
	rm -rf *.o *~ \#*\#

clean: tidy
