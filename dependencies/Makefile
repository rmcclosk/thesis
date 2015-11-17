GNU = ftp://gnu.mirror.iweb.com
WD = $(shell pwd)

all: src \
     bin/tar \
     lib/liblzma.so \
     bin/automake \
     bin/autoconf \
     bin/libtool \
     lib/libcurl.so \
     bin/git \
     bin/pkg-config \
     bin/flex \
     bin/bison \
     lib/libgsl.so \
     lib/libigraph.so \
     lib/libJudy.so \
     lib/libyaml.so \
     lib/libcheck.so

src:
	mkdir $@

bin/tar: src/tar-1.28.tar.gz
lib/liblzma.so: src/xz-5.2.2.tar.gz
lib/libcurl.so: src/curl-7.45.0.tar.gz
bin/automake: src/automake-1.15.tar.xz
bin/autoconf: src/autoconf-2.69.tar.xz
bin/libtool: src/libtool-2.4.6.tar.xz
bin/flex: src/flex-2.6.0.tar.bz2
bin/bison: src/bison-3.0.4.tar.gz
lib/libgsl.so: src/gsl-2.1.tar.gz
lib/libJudy.so: src/judy-1.0.5.tar.gz
lib/libyaml.so: src/yaml-0.1.5.tar.gz
lib/libcheck.so: src/check-0.10.0.tar.gz

bin/pkg-config: src/pkg-config-0.29.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure --prefix=$(WD) --with-internal-glib && make && make install

lib/libigraph.so: src/igraph
	cd $^ && ./bootstrap.sh && ./configure --prefix=$(WD) --enable-tls && make && make install

bin/git: src/git-2.6.3.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && make prefix=$(WD) all && make prefix=$(WD) install

bin/% lib/%:
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure --prefix=$(WD) && make && make install

src/yaml-0.1.5.tar.gz:
	curl -z $@ -o $@ http://pyyaml.org/download/libyaml/yaml-0.1.5.tar.gz

src/igraph:
	cd $(dir $@) && git clone https://github.com/igraph/igraph.git

src/judy-1.0.5.tar.gz:
	curl -z $@ -o $@ http://iweb.dl.sourceforge.net/project/judy/judy/Judy-1.0.5/Judy-1.0.5.tar.gz

src/git-2.6.3.tar.gz:
	curl -z $@ -o $@ https://codeload.github.com/git/git/tar.gz/v2.6.3

src/pkg-config-0.29.tar.gz:
	curl -z $@ -o $@ http://pkgconfig.freedesktop.org/releases/pkg-config-0.29.tar.gz

src/check-0.10.0.tar.gz:
	curl -z $@ -o $@ http://skylineservers.dl.sourceforge.net/project/check/check/0.10.0/check-0.10.0.tar.gz

src/xz-5.2.2.tar.gz:
	curl -z $@ -o $@ http://tukaani.org/xz/xz-5.2.2.tar.gz

src/curl-7.45.0.tar.gz:
	curl -z $@ -o $@ http://curl.haxx.se/download/curl-7.45.0.tar.gz

src/flex-2.6.0.tar.bz2:
	curl -z $@ -o $@ http://iweb.dl.sourceforge.net/project/flex/flex-2.6.0.tar.bz2

src/%.tar.xz src/%.tar.gz:
	curl -z $@ -o $@ --disable-epsv $(GNU)/$(word 1, $(subst -, , $(notdir $@)))/$(notdir $@)
