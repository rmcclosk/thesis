GNU = ftp://gnu.mirror.iweb.com
CRAN = http://cran.stat.sfu.ca
WD = $(shell pwd)

all: dirs \
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
     bin/R \
     lib/libgsl.so \
     lib/libigraph.so \
     lib/libJudy.so \
     lib/libyaml.so \
     lib/libcheck.so \
     lib/libxml2.so \
	 bin/python3 \
     lib/python3/site-packages/yaml \
     R-pkgs

dirs:
	mkdir -p src R/library

R-pkgs: R/library/devtools \
        R/library/kernlab \
        R/library/ggplot2 \
        R/library/yaml \
        R/library/ape \
        R/library/igraph \
        R/library/RColorBrewer \
        R/library/data.table \
        R/library/reshape2 \
        R/library/VGAM

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
bin/python3: src/Python-3.5.0.tar.xz
lib/python3/site-packages/yaml: src/PyYAML-3.11.tar.gz

bin/pkg-config: src/pkg-config-0.29.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure --prefix=$(WD) --with-internal-glib && make && make install

lib/libigraph.so: src/igraph
	cd $^ && ./bootstrap.sh && ./configure --prefix=$(WD) --enable-tls && make && make install

bin/git: src/git-2.6.3.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && make prefix=$(WD) all && make prefix=$(WD) install

bin/R: src/R-3.2.2.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure --prefix=$(WD) --enable-R-shlib && make && make install

lib/libxml2.so: src/libxml2-2.9.2.tar.gz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure --prefix=$(WD) --without-python && make && make install

lib/python3/site-packages/%:
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && python3 setup.py install

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

src/R-3.2.2.tar.gz:
	curl -z $@ -o $@ $(CRAN)/src/base/R-3/R-3.2.2.tar.gz

src/libxml2-2.9.2.tar.gz:
	curl -z $@ -o $@ --disable-epsv ftp://xmlsoft.org/libxml2/libxml2-2.9.2.tar.gz

src/Python-3.5.0.tar.xz:
	curl -z $@ -o $@ https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tar.xz

src/PyYAML-3.11.tar.gz:
	curl -z $@ -o $@ http://pyyaml.org/download/pyyaml/PyYAML-3.11.tar.gz

src/%.tar.xz src/%.tar.gz:
	curl -z $@ -o $@ --disable-epsv $(GNU)/$(word 1, $(subst -, , $(notdir $@)))/$(notdir $@)

R/library/%:
	R -e 'options(repos=structure(c(CRAN="$(CRAN)")))' -e 'install.packages("$*")'
