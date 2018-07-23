# use ubuntu image as base
from ubuntu:17.10

# set the working directory to DualSPHysics
WORKDIR ./DualSPHysics

# copy the current directory in this folder
ADD ./bin/linux/ /DualSPHysics

# install build essentials
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils libgomp1

# make files executable
RUN chmod +x ./BoundaryVTK4_linux64
RUN chmod +x ./ComputeForces4_linux64
RUN chmod +x ./DualSPHysics4.2CPU_linux64
RUN chmod +x ./DualSPHysics4.2_linux64
RUN chmod +x ./FloatingInfo4_linux64
RUN chmod +x ./GenCase4_MkWord_linux64
RUN chmod +x ./GenCase4_linux64
RUN chmod +x ./IsoSurface4_linux64
RUN chmod +x ./MeasureBoxes4_linux64
RUN chmod +x ./MeasureTool4_linux64
RUN chmod +x ./PartVTK4_linux64
RUN chmod +x ./PartVTKOut4_linux64
