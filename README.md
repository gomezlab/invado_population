invado_population
=================

Quantification of Invadopodia activity in cell populations

This repository contains a set of MATLAB scripts designed to process images taken of metastatic cancer cells degrading labeled ECM. After finding and tracking individual cells, the system quantifies the amount of degradation underneath each cell over time. The system has automated systems for flat-field and photobleach correcting the labeled ECM images.

Data Organization and Processing
================================

The system expects your images to be organized in a specific manner. Basicially, each field of view has a folder, with subfolders for the ECM and cell marker images. Each of these field of view folders also must contain a .cfg file which points to the master config file that tells the computer where to find the rest of the image sets.
