FROM laristra/flecsi:fedora_serial

ARG MINIMAL
ARG COVERAGE
ARG CC
ARG CXX
ARG CXXFLAGS
ARG BUILD_TYPE
ARG INTERNAL_FLECSI

#for coverage
ENV CI true
ENV TRAVIS true
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST 
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT

RUN rm -rf /home/flecsi/.ccache
COPY flecsale /home/flecsi/flecsale
COPY ccache/ /home/flecsi/.ccache
USER root
RUN chown -R flecsi:flecsi /home/flecsi/flecsale /home/flecsi/.ccache
USER flecsi

WORKDIR /home/flecsi/flecsale
RUN mkdir build

WORKDIR build
RUN ccache -z
RUN cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          ${MINIMAL:+-DCMAKE_DISABLE_FIND_PACKAGE_EXODUSII=ON \
          -DCMAKE_DISABLE_FIND_PACKAGE_SCOTCH=ON \
          -DCMAKE_DISABLE_FIND_PACKAGE_METIS=ON  \
          -DCMAKE_DISABLE_FIND_PACKAGE_LAPACKE=ON} \
          -DCMAKE_DISABLE_FIND_PACKAGE_FleCSI=${INTERNAL_FLECSI} \
          -DENABLE_UNIT_TESTS=ON \
          -DENABLE_DOXYGEN=ON \
          ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} ..
RUN make VERBOSE=1 -j2
RUN ccache -s
RUN make test ARGS="-V"
RUN make install DESTDIR=${PWD}
RUN make doxygen flecsi.doxygen
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov; \
  fi; \
fi && cd -
USER root
RUN make install
USER flecsi
WORKDIR /home/flecsi
