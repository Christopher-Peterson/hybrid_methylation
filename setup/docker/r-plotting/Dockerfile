FROM crpeters/r-tidyverse-optparse:4.2.1

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y \
 libx11-dev libcairo2-dev libxt-dev && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*


RUN install2.r -e ggforce ggExtra ggh4x \
  ggpointdensity GGally geomtextpath \
  ggrepel ggdist ggtext patchwork \
  cowplot ggalt gganimate ggraph \
  ggbump ggalluvial ggridges \
  gggenes ComplexUpset ggpattern \
  ggupset ggside ggTimeSeries \
  gghalves ggdensity ggcorrplot \
  ggmuller directlabels ggheatmap \
  ggimg ggarchery Cairo ggrastr ragg