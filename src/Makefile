all: libfoo.so libbar.so libharmonic_mean.so

libharmonic_mean.so: harmonic_mean.c
	$(CC) $(CFLAGS) -lm -shared -fPIC -o libharmonic_mean.so harmonic_mean.c

libfoo.so: foo.c
	$(CC) $(CFLAGS) -shared -fPIC -o libfoo.so foo.c

libbar.so: bar.c
	$(CC) $(CFLAGS) -shared -fPIC -o libbar.so bar.c

tidy:
	rm -rf *.o *~ \#*\#

clean: tidy
