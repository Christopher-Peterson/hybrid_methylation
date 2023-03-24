FROM niemasd/samtools:1.16.1

RUN wget https://github.com/FelixKrueger/SNPsplit/archive/refs/tags/0.5.0.tar.gz && \
	tar xvzf 0.5.0.tar.gz && \
	rm 0.5.0.tar.gz && \
	cd SNPsplit-0.5.0 && \
	chmod a+x tag2sort SNPsplit_genome_preparation SNPsplit

ENV PATH=$PATH:/usr/src/SNPsplit-0.5.0/