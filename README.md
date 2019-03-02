# cracking_codes

Back up and examples of various codes written to solve cracking problems

cracking_uniform contains the code for cracking in uniform films
cracking_drying contains coupled codes to solve for fracture in drying droplets

Codes solve the governing equations using finite elements, within the FeNiCs enviroment

nf_ denotes codes written for NextFlow pipeline, where solutions are read using YAML

To create singularity shell, run the following command as adminstrator
 sudo singularity build fen_crack.def fen_crack.sif

Relevant Docker image can be found within the fen_crack.def
