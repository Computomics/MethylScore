Bootstrap:docker
From:conda/miniconda2

%files
    environment.yml /

%post
    /usr/local/bin/conda env update -n root -f /environment.yml
    /usr/local/bin/conda clean -a
    mkdir -p /lustre
