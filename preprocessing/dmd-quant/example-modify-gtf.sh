

#!/usr/bin/env bash

OTHER="\
ENSMUST00000113992|\
ENSMUST00000123308|\
ENSMUST00000113991|\
ENSMUST00000127295|\
ENSMUST00000128983|\
ENSMUST00000132333|\
ENSMUST00000139998|\
ENSMUST00000141261|\
ENSMUST00000146331|\
ENSMUST00000147740|\
ENSMUST00000149433|\
ENSMUST00000156107\
"

zcat Mus_musculus.GRCm39.104.gtf.gz | grep -Pv "$OTHER" > Mus_musculus_dmd_2transcripts140.gtf
