FROM bids/base_fsl

MAINTAINER Cameron Craddock <cameron.craddock@childmind.org>

RUN curl -sSLO http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    bash Miniconda2-latest-Linux-x86_64.sh -b && \
    /root/miniconda2/bin/conda update -yq conda && \
    /root/miniconda2/bin/conda install -y -c conda-forge scipy nipype pandas patsy && \
    /root/miniconda2/bin/conda clean -y --all


# add all of the environment variables for fsl

ENV FSLDIR /usr/share/fsl/5.0
ENV FSLOUTPUTTYPE NIFTI_GZ
ENV FSLMULTIFILEQUIT TRUE
ENV FSLTCLSH /usr/bin/tclsh
ENV FSLWISH /usr/bin/wish
ENV FSLBROWSER /etc/alternatives/x-www-browser
ENV LD_LIBRARY_PATH /usr/lib/fsl/5.0:${LD_LIBRARY_PATH}
ENV PATH ${FSLDIR}/bin:/root/miniconda2/bin:${PATH}

COPY create_flame_model_files.py /code/create_flame_model_files.py
COPY run_new_2.py /code/run_new_2.py

ENTRYPOINT ["/code/run_new_2.py"]
