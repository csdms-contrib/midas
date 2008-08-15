
GFORTRAN=gfortran
FFLAGS=
PROG=midas
VERSION=0.1
SOURCES=fldta.f entrain.f entrainh.f settle.f svela.f turb.f bedload.f logdist.f susp.f midas.f \
        erodep.f drange.f dfifty.f drangeh.f susp1.f elwidcalc.f wrfldta.f writerstrt.f table.f \
        toprint.f result.f plotfiles.f discr.f readrstrt.f gsdist.f
OBJS=${SOURCES:.f=.o}

all: ${PROG}

${PROG}: ${SOURCES}
	${GFORTRAN} ${FFLAGS} -o ${PROG} ${SOURCES}

dist:
	@mkdir ${PROG}-${VERSION}
	@cp ${SOURCES} ${PROG}-${VERSION}
	@cp Makefile ${PROG}-${VERSION}
	@tar cvfz ${PROG}-${VERSION}.tar.gz ${PROG}-${VERSION} 
	@rm -rf ${PROG}-${VERSION}

clean:
	@rm -f ${PROG} ${OBJS} core

