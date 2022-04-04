---
title: "Matlab Code Documentation"
permalink: /code-doc/
---


# Design of code structure base on functionality

This is where I can think about using Unified modeling language (UML).

1. feature identification
  - tool: regionprops
    - general requirements: threshold for feature detection, 
    - data specific requirements: Low level satellite data would require info. about cloud mask and how to deal with cloud mask during the feature identification; high level gap-filled data would not need this info.

  - data storing structure: in a matlab structure that can be thought as a huge dictionary in python. 

  - utility functions: 
    - take information out of the stored data units
    - aggregate information needed in a better way (e.g., blob characteristics)
    - find spatial probability distribution of the features and view them in a standard way.
    - ?
  

2. extract cloudiness over features
  - utility functions:
    - cloudiness computation from raw cloud masks (different ways to compute cloudiness)
    - coordinate transformation (yaxis align with the wind direction)
      - require wind information (allow flexibility to take in wind information from different height levels)
      - coordinate is feature centered and normalized based on the size of the features
        - provide two different options for normalization: ellipse versus cirle

3. explore parameter dependence
   - utility functions:
     - sorting data (features & cloudiness) by different criteria
       - 1. Atmospheric conditions (surface wind and lower tropospheric stability)
       - 2. feature characteristics (size, strength, etc)
       - 2. surface wind
       - EOF to find the leading patterns that are associated with the main variance.

# Feature detection:
|-----------------------------|
| script name | Purpose       |
| ----------- | ------------- |
| 


# probability analysis:
|-----------------------------|
| script name | Purpose       |
| ----------- | ------------- |





