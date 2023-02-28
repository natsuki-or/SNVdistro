#!/usr/bin/env python3.6

import getuniprot as gu

id = ['003R_FRG3G','C0M2l1_BPMS2','A0A7C1M9G7_9ARCH','005L_IIV3','CDO1_HUMAN','zzz','006L_IIV6','P21816']
output = gu.getuniprot(id)
for line in output:
    print (line, end="")
