#!/usr/bin/env python3
import re
dict1={}
with open('all_tail_changeID_1207.protein.fa_70.clstr','r') as f1:
    for line in f1.readlines():
        line=line.strip()
        if re.match('>',line):
            cluster=line.split('>')[-1]
           # dict1[ID]=[]
        else:
            phage=line.split(', >')[-1].split('...')[0]
            dict1[phage]=cluster
#list1=['PHKP9582-1_YP_009191451.1','PHKP9828_A1_YP_009201600.1','PH-KP8949-A_YP_004010030.1','PH-KP9115-1_ALT58490.1','PH-KP9202-B1_ALT58490.1','PH-KP9247-C1_AZF88836.1','PH-KP9072-1_ALT58490.1','PH-KP8870-A1_ALT58490.1','PHKP4049_AZF88836.1','PH-KP7215-A_AWN07165.1','PH-KP8565_AOZ65513.1','PH-KP7569_ALT58490.1','PHKP9088-A_ALT58490.1','PH-KP7168_pdb|5MU4|D','PH-KP5137_ALT58490.1','PH-KP5556_AWN07207.1','PHKP9088-A_QCG76439.1','PHKP9495-AB1_YP_009146416.1','PH-KP7524_ASW27611.1','PH-KP9841-A1_ANN86331.1','PHKP9822_B1_YP_009626305.1','PH-KP9353-A1_ANN86331.1','PH-KP9091-1_ANN86331.1','PH-KP9887-1_ASW27611.1','PHKP9346_A1_YP_009604496.1','PHKP9778_C1_YP_009609870.1','PH-KP5996_AUG87747.1','PH-KP7526_AUG87747.1','PH-KP5295_ATS92567.1','PH-KP7654_YP_009615313.1','PHKP9025-A1_YP_007392887.1','PH-KP7338_AYP28214.1','PH-KP7148_AYP28214.1','PH-KP5431_AYP28214.1','PHKP9060-1_AVI03135.1','PH-KP9083-A_YP_009284925.1','PH-KP9595-A1_YP_009284925.1','PHKP9721-A_AXQ68098.1','PH-KP4547_YP_009195382.1','PHKP9773_C1_ARM70345.1','PH-KP9140-1_AXQ68098.1','PH-KP9332_A1_YP_009226010.1','PH-KP9716-A_YP_009226010.1','PH-KP9349_A1_QBQ72942.1','PHKP9495-AB1_BAW85700.1','PH-KP5130-A_ALT58497.1','PH-KP9438_A1_QDB71175.1','PH-KP5493-A_AUV57596.1','PH-KP6604_AUV57596.1','PH-KP5493-A_YP_009288925.1','PH-KP9841-A1_YP_009626439.1','PHKP9822_B1_YP_007348839.1','PH-KP9091-1_YP_009194451.1','PH-KP6604_AWY07055.1','PH-KP5493-A_AWY07055.1','PH-KP5119_BAW85701.1','PH-KP9438-B1_BAW85701.1','PH-KP8956-AB1_BAW85701.1','PH-KP7102_BAW85701.1','PHKP9283_B1_AWN07078.1','PHKP4049_AWN07078.1','PH-KP9247-C1_AWN07078.1','PH-KP8956-AB1_YP_009153187.1','PH-KP5119_YP_009153187.1','PH-KP7102_YP_009153187.1','PH-KP4829_YP_009153187.1','PH-KP5066_AWD90221.1','PH-KP5619_YP_009199639.1','PH-KP9850-B1_YP_007348840.1','PH-KP9841-A1_YP_007348840.1','PH-KP9091-1_YP_003580044.1','PH-KP9353-A1_YP_009607401.1','PHKP9495-AB1_YP_009146411.1','PH-KP9380-A1_AUG87917.1','PH-KP6031_AUG87917.1','PH-KP6221_AUG87917.1','PH-KP5996_YP_008532017.1','PHKP9261-A1_YP_008532017.1','PH-KP7138_YP_008532017.1','PH-KP5619_YP_009199654.1','PH-KP9310-A1_QBQ72941.1','PHKP9434_B1_ASW27652.1','PH-KP5150_AXC43025.1','PHKP9947-AB1_AWP45426.1','PHKP9429_B1_AXC43025.1','PHKP9947-B1_AXC43025.1','PH-KP5130-A_QBQ72011.1','PHKP-9534-A1_YP_009302749.1','PHKP9310_B1_hypothesis_protein_31876_32043','PH-KP6604_AUV57507.1','PH-KP9074-1_ALJ98202.1','PH-KP4956-A_ALJ98202.1','PH-KP8949-A_AVR55358.1','PH-KP5287_AUG87959.1','PH-KP7215-A_ASZ78307.1','PH-KP5556_APZ82804.1','PH-KP9202-B1_YP_009204835.1','PH-KP5137_AOT28172.1','PH-KP9115-1_AWN07083.1','PH-KP9072-1_ALT58497.1','PH-KP8870-A1_ALT58497.1','PH-KP9233-A1_ALT58497.1','PH-KP7168_YP_002003830.1','PH-KP7569_ALT58497.1','PHKP9450-A1_AZF88843.1','PHKP9450-B1_AZF88843.1','PH-KP9247-C1_AZF88843.1','PH-KP8565_AZF88843.1','PH-KP9074-1_ALJ98203.1','PHKP9947-AS1_SPF82155.1','PH-KP9628-B1_AUV62702.1','PHKP9721-A_ASW27458.1','PH-KP9140-1_AXQ68100.1','PH-KP9298-A1_AXQ68100.1','PHKP9060-1_AUG87959.1','PH-KP6208_YP_008532048.1','PH-KP7148_AYP28213.1','PH-KP9332_A1_AYP28213.1','PHKP9716_A1_AYP28213.1','PH-KP9716-A_AYP28213.1','PH-KP7639_AUG87748.1','PH-KP9595-A1_BAW85698.1','PH-KP9083-A_BAW85698.1','PH-KP7338_AYP28213.1','PHKP9947-AB1_QAX91947.1','PHKP9947-B1_QAX91947.1','PHKP9429_B1_QAX91947.1','PHKP8949_A1_YP_004010112.1','PH-KP8949-A_YP_004010112.1','PHKP9822_B1_YP_009626500.1','PH-KP9091-1_YP_007348899.1','PH-KP9353-A1_YP_007348899.1','PHKP9298-B1_YP_007348899.1','PHKP9748-A1_ASW27611.1','PH-KP7524_ASW27611.1','PH-KP9628-B1_ASW27611.1','PH-KP5087_ASW27611.1','PH-KP9091-1_YP_009626502.1','PH-KP9353-A1_YP_007348901.1','PH-KP9841-A1_YP_009626502.1','PHKP9822_B1_YP_009626502.1','PH-KP7760_YP_009188359.1','PH-KP5048_YP_009098379.1','PH-KP9566-1_YP_009188359.1','PH-KP8956-AB1_AUE23419.1','PH-KP7102_BAW85695.1','PH-KP5493-A_AWY07059.1','PH-KP6604_AWY07059.1','PH-KP5493-A_AWY07149.1','PH-KP5493-A_YP_009288927.1','PH-KP5431_YP_008532048.1','PHKP9429_B1_QBJ02991.1','PHKP9495-AB1_YP_009146412.1','PH-KP5087_AUV62702.1','PH-KP7102_BAW85692.1','PHKP8949_A1_YP_004010021.1','PHKP9429_C1_YP_005098088.1','PHKP9947-AB1_YP_005098088.1','PH-KP7569_BBF66885.1','PH-KP7168_YP_002003825.1','PH-KP8565_BBF66862.1','PH-KP7215-A_AWN07120.1','PH-KP5137_AWN07078.1','PHKP9088-A_AWN07078.1','PH-KP9115-1_AWN07078.1','PH-KP9202-B1_QCG76444.1','PH-KP9233-A1_YP_003347550.1','PH-KP8870-A1_ARB12447.1','PH-KP5556_AOT28167.1','PH-KP6604_AWY07147.1','PH-KP5493-A_AUV57597.1','PH-KP9384_A1_ASW27497.1','PH-KP6221_AUG87747.1','PH-KP9380-A1_AUG87747.1','PHKP9242_A1_AUG87747.1','PH-KP6031_AUG87747.1','PHKP-9529-1_AUG87747.1','PH-KP9269-A1_AUG87747.1','PH-KP7526_AUG87747.1','PHKP9467-1_AUG87747.1','PH-KP9052-A1_AUG87747.1','PH-KP5996_AUG87747.1','PHKP9261-A1_AUG87747.1','PH-KP7138_AUG87747.1','PHKP9298-B1_YP_009607461.1','PH-KP9353-A1_YP_009194513.1','PH-KP9091-1_YP_009194513.1','PH-KP9841-A1_YP_009194513.1','PHKP9822_B1_YP_009626501.1','PH-KP9091-1_YP_009194450.1','PH-KP9353-A1_YP_009194450.1','PH-KP9841-A1_YP_009626438.1','PHKP9822_B1_YP_007348838.1','PH-KP9091-1_YP_009607397.1','PHKP9822_B1_YP_009607397.1','PH-KP9841-A1_YP_007348836.1','PH-KP5295_ATS92566.1','PH-KP7654_YP_009615317.1','PHKP9310_B1_hypothesis_protein_31459_31683','PHKP9261-A1_ATS92567.1','PH-KP6031_AYP69336.1','PH-KP7138_AUG87749.1','PH-KP6221_YP_007007682.1','PH-KP9380-A1_AUG87749.1','PH-KP5066_BAT32505.1','PH-KP9091-1_YP_009607107.1','PH-KP9841-A1_YP_009607107.1','PH-KP9353-A1_YP_009607107.1','PHKP9822_B1_YP_009607107.1','PH-KP9353-A1_YP_007348897.1','PHKP9298-B1_YP_007348897.1','PH-KP9091-1_YP_007348897.1','PH-KP9850-B1_YP_007348897.1','PH-KP9841-A1_YP_007348897.1','PHKP9822_B1_YP_003580107.1','PHKP9593-1_pdb|4UW8|I','PH-KP9247-C1_AZF88843.1','PH-KP5493-A_YP_009288844.1','PHKP9298-B1_YP_007348898.1','PH-KP9841-A1_YP_007348898.1','PHKP9822_B1_YP_003580108.1','PH-KP9353-A1_YP_007348898.1','PH-KP9091-1_YP_007348898.1','PHKP9582-1_QBZ71284.1','PH-KP6017_YP_007008117.1','PH-KP9438_A1_YP_003347555.1','PH-KP5130-A_YP_003347555.1','PH-KP7102_YP_007007689.1','PH-KP5119_YP_007007689.1','PH-KP5119_YP_007007686.1','PH-KP7102_ATS92567.1','PH-KP8956-AB1_YP_007007686.1','PH-KP9072-1_AOZ65570.1','PH-KP9233-A1_AOZ65570.1','PH-KP7569_AOZ65570.1','PH-KP9115-1_AOT28173.1','PH-KP5137_AWN07084.1','PH-KP9310-A1_ACH72967.1','PH-KP7673_AUG87749.1','PH-KP9385-A1_QDB71140.1','PH-KP9438_A1_YP_009609192.1','PHKP9773_C1_QBQ72011.1','PH-KP9072-1_YP_009190990.1','PH-KP7338_AXN53795.1','PH-KP9332_A1_AXN53795.1','PH-KP9716-A_AXN53795.1','PH-KP5431_AXN53795.1','PH-KP9140-1_AXN53795.1','PH-KP7148_AXQ67940.1','PH-KP4547_YP_009284929.1','PH-KP9595-A1_YP_009284929.1','PHKP9773_C1_YP_009284929.1','PH-KP7639_YP_009284929.1','PHKP9721-A_YP_009284929.1','PH-KP9083-A_YP_009284929.1','PH-KP9064-1_YP_009188359.1','PHKP9471-A1_YP_009199929.1','PH-KP5406_YP_009199929.1','PHKP9346_A1_YP_009609871.1','PHKP9778_C1_YP_009592223.1','PH-KP9072-1_hypothesis_protein_1000_1458','PH-KP9332_A1_tr|A0A286MN37|A0A286MN37_9CAUD','PHKP9822_B1_YP_007348840.1','PH-KP6017_AUX83696.1','PHKP4049_hypothesis_protein_39750_39962','PHKP9283_B1_AZF88843.1','PHKP4049_AZF88843.1']
list1=['PHKP9582-1_YP_009191451.1','PHKP9828_A1_YP_009201600.1','PH-KP8949-A_YP_004010030.1','PH-KP9115-1_ALT58490.1','PH-KP9202-B1_ALT58490.1','PH-KP9247-C1_AZF88836.1','PH-KP9072-1_ALT58490.1','PH-KP8870-A1_ALT58490.1','PHKP4049_AZF88836.1','PH-KP7215-A_AWN07165.1','PH-KP8565_AOZ65513.1','PH-KP7569_ALT58490.1','PHKP9088-A_ALT58490.1','PH-KP7168_pdb|5MU4|D','PH-KP5137_ALT58490.1','PH-KP5556_AWN07207.1','PHKP9088-A_QCG76439.1','PHKP9495-AB1_YP_009146416.1','PH-KP7524_S_ASW27611.1','PH-KP9841-A1_ANN86331.1','PHKP9822_B1_YP_009626305.1','PH-KP9353-A1_ANN86331.1','PH-KP9091-1_ANN86331.1','PH-KP9887-1_ASW27611.1','PHKP9346_A1_YP_009604496.1','PHKP9778_C1_YP_009609870.1','PH-KP5996_s_AUG87747.1','PH-KP7526_s_AUG87747.1','PH-KP5295_ATS92567.1','PH-KP7654_YP_009615313.1','PHKP9025-A1_YP_007392887.1','PH-KP7338_AYP28214.1','PH-KP7148_AYP28214.1','PH-KP5431_AYP28214.1','PHKP9060-1_AVI03135.1','PH-KP9083-A_YP_009284925.1','PH-KP9595-A1_YP_009284925.1','PHKP9721-A_AXQ68098.1','PH-KP4547_YP_009195382.1','PHKP9773_C1_ARM70345.1','PH-KP9140-1_AXQ68098.1','PH-KP9332_A1_YP_009226010.1','PH-KP9716-A_YP_009226010.1','PH-KP9349_A1_QBQ72942.1','PHKP9495-AB1_BAW85700.1','PH-KP5130-A_ALT58497.1','PH-KP9438_A1_QDB71175.1','PH-KP5493-A_AUV57596.1','PH-KP6604_AUV57596.1','PH-KP5493-A_YP_009288925.1','PH-KP9841-A1_YP_009626439.1','PHKP9822_B1_YP_007348839.1','PH-KP9091-1_YP_009194451.1','PH-KP6604_AWY07055.1','PH-KP5493-A_AWY07055.1','PH-KP5119_BAW85701.1','PH-KP9438-B1_BAW85701.1','PH-KP8956-AB1_BAW85701.1','PH-KP7102_BAW85701.1','PHKP9283_B1_AWN07078.1','PHKP4049_AWN07078.1','PH-KP9247-C1_AWN07078.1','PH-KP8956-AB1_YP_009153187.1','PH-KP5119_YP_009153187.1','PH-KP7102_YP_009153187.1','PH-KP4829_YP_009153187.1','PH-KP5066_AWD90221.1','PH-KP5619_YP_009199639.1','PH-KP9850-B1_YP_007348840.1','PH-KP9841-A1_YP_007348840.1','PH-KP9091-1_YP_003580044.1','PH-KP9353-A1_YP_009607401.1','PHKP9495-AB1_YP_009146411.1','PH-KP9380-A1_AUG87917.1','PH-KP6031_AUG87917.1','PH-KP6221_AUG87917.1','PH-KP5996_YP_008532017.1','PHKP9261-A1_YP_008532017.1','PH-KP7138_YP_008532017.1','PH-KP5619_YP_009199654.1','PH-KP9310-A1_QBQ72941.1','PHKP9434_B1_ASW27652.1','PH-KP5150_AXC43025.1','PHKP9947-AB1_AWP45426.1','PHKP9429_B1_AXC43025.1','PHKP9947-B1_AXC43025.1','PH-KP5130-A_QBQ72011.1','PHKP-9534-A1_YP_009302749.1','PHKP9310_B1_hypothesis_protein_31876_32043','PH-KP6604_AUV57507.1','PH-KP9074-1_ALJ98202.1','PH-KP4956-A_ALJ98202.1','PH-KP8949-A_AVR55358.1','PH-KP5287_AUG87959.1','PH-KP7215-A_ASZ78307.1','PH-KP5556_APZ82804.1','PH-KP9202-B1_YP_009204835.1','PH-KP5137_AOT28172.1','PH-KP9115-1_AWN07083.1','PH-KP9072-1_ALT58497.1','PH-KP8870-A1_ALT58497.1','PH-KP9233-A1_ALT58497.1','PH-KP7168_YP_002003830.1','PH-KP7569_ALT58497.1','PHKP9450-A1_AZF88843.1','PHKP9450-B1_AZF88843.1','PH-KP9247-C1_L_AZF88843.1','PH-KP8565_AZF88843.1','PH-KP9074-1_ALJ98203.1','PHKP9947-AS1_SPF82155.1','PH-KP9628-B1_AUV62702.1','PHKP9721-A_ASW27458.1','PH-KP9140-1_AXQ68100.1','PH-KP9298-A1_AXQ68100.1','PHKP9060-1_AUG87959.1','PH-KP6208_YP_008532048.1','PH-KP7148_AYP28213.1','PH-KP9332_A1_AYP28213.1','PHKP9716_A1_AYP28213.1','PH-KP9716-A_AYP28213.1','PH-KP7639_AUG87748.1','PH-KP9595-A1_BAW85698.1','PH-KP9083-A_BAW85698.1','PH-KP7338_AYP28213.1','PHKP9947-AB1_QAX91947.1','PHKP9947-B1_QAX91947.1','PHKP9429_B1_QAX91947.1','PHKP8949_A1_YP_004010112.1','PH-KP8949-A_YP_004010112.1','PHKP9822_B1_YP_009626500.1','PH-KP9091-1_YP_007348899.1','PH-KP9353-A1_YP_007348899.1','PHKP9298-B1_YP_007348899.1','PHKP9748-A1_ASW27611.1','PH-KP7524_L_ASW27611.1','PH-KP9628-B1_ASW27611.1','PH-KP5087_ASW27611.1','PH-KP9091-1_YP_009626502.1','PH-KP9353-A1_YP_007348901.1','PH-KP9841-A1_YP_009626502.1','PHKP9822_B1_YP_009626502.1','PH-KP7760_YP_009188359.1','PH-KP5048_YP_009098379.1','PH-KP9566-1_YP_009188359.1','PH-KP8956-AB1_AUE23419.1','PH-KP7102_BAW85695.1','PH-KP5493-A_AWY07059.1','PH-KP6604_AWY07059.1','PH-KP5493-A_AWY07149.1','PH-KP5493-A_YP_009288927.1','PH-KP5431_YP_008532048.1','PHKP9429_B1_QBJ02991.1','PHKP9495-AB1_YP_009146412.1','PH-KP5087_AUV62702.1','PH-KP7102_BAW85692.1','PHKP8949_A1_YP_004010021.1','PHKP9429_C1_YP_005098088.1','PHKP9947-AB1_YP_005098088.1','PH-KP7569_BBF66885.1','PH-KP7168_YP_002003825.1','PH-KP8565_BBF66862.1','PH-KP7215-A_AWN07120.1','PH-KP5137_AWN07078.1','PHKP9088-A_AWN07078.1','PH-KP9115-1_AWN07078.1','PH-KP9202-B1_QCG76444.1','PH-KP9233-A1_YP_003347550.1','PH-KP8870-A1_ARB12447.1','PH-KP5556_AOT28167.1','PH-KP6604_AWY07147.1','PH-KP5493-A_AUV57597.1','PH-KP9384_A1_ASW27497.1','PH-KP6221_AUG87747.1','PH-KP9380-A1_AUG87747.1','PHKP9242_A1_AUG87747.1','PH-KP6031_AUG87747.1','PHKP-9529-1_AUG87747.1','PH-KP9269-A1_AUG87747.1','PH-KP7526_L_AUG87747.1','PHKP9467-1_AUG87747.1','PH-KP9052-A1_AUG87747.1','PH-KP5996_L_AUG87747.1','PHKP9261-A1_AUG87747.1','PH-KP7138_AUG87747.1','PHKP9298-B1_YP_009607461.1','PH-KP9353-A1_YP_009194513.1','PH-KP9091-1_YP_009194513.1','PH-KP9841-A1_YP_009194513.1','PHKP9822_B1_YP_009626501.1','PH-KP9091-1_YP_009194450.1','PH-KP9353-A1_YP_009194450.1','PH-KP9841-A1_YP_009626438.1','PHKP9822_B1_YP_007348838.1','PH-KP9091-1_YP_009607397.1','PHKP9822_B1_YP_009607397.1','PH-KP9841-A1_YP_007348836.1','PH-KP5295_ATS92566.1','PH-KP7654_YP_009615317.1','PHKP9310_B1_hypothesis_protein_31459_31683','PHKP9261-A1_ATS92567.1','PH-KP6031_AYP69336.1','PH-KP7138_AUG87749.1','PH-KP6221_YP_007007682.1','PH-KP9380-A1_AUG87749.1','PH-KP5066_BAT32505.1','PH-KP9091-1_YP_009607107.1','PH-KP9841-A1_YP_009607107.1','PH-KP9353-A1_YP_009607107.1','PHKP9822_B1_YP_009607107.1','PH-KP9353-A1_YP_007348897.1','PHKP9298-B1_YP_007348897.1','PH-KP9091-1_YP_007348897.1','PH-KP9850-B1_YP_007348897.1','PH-KP9841-A1_YP_007348897.1','PHKP9822_B1_YP_003580107.1','PHKP9593-1_pdb|4UW8|I','PH-KP9247-C1_s_AZF88843.1','PH-KP5493-A_YP_009288844.1','PHKP9298-B1_YP_007348898.1','PH-KP9841-A1_YP_007348898.1','PHKP9822_B1_YP_003580108.1','PH-KP9353-A1_YP_007348898.1','PH-KP9091-1_YP_007348898.1','PHKP9582-1_QBZ71284.1','PH-KP6017_YP_007008117.1','PH-KP9438_A1_YP_003347555.1','PH-KP5130-A_YP_003347555.1','PH-KP7102_YP_007007689.1','PH-KP5119_YP_007007689.1','PH-KP5119_YP_007007686.1','PH-KP7102_ATS92567.1','PH-KP8956-AB1_YP_007007686.1','PH-KP9072-1_AOZ65570.1','PH-KP9233-A1_AOZ65570.1','PH-KP7569_AOZ65570.1','PH-KP9115-1_AOT28173.1','PH-KP5137_AWN07084.1','PH-KP9310-A1_ACH72967.1','PH-KP7673_AUG87749.1','PH-KP9385-A1_QDB71140.1','PH-KP9438_A1_YP_009609192.1','PHKP9773_C1_QBQ72011.1','PH-KP9072-1_YP_009190990.1','PH-KP7338_AXN53795.1','PH-KP9332_A1_AXN53795.1','PH-KP9716-A_AXN53795.1','PH-KP5431_AXN53795.1','PH-KP9140-1_AXN53795.1','PH-KP7148_AXQ67940.1','PH-KP4547_YP_009284929.1','PH-KP9595-A1_YP_009284929.1','PHKP9773_C1_YP_009284929.1','PH-KP7639_YP_009284929.1','PHKP9721-A_YP_009284929.1','PH-KP9083-A_YP_009284929.1','PH-KP9064-1_YP_009188359.1','PHKP9471-A1_YP_009199929.1','PH-KP5406_YP_009199929.1','PHKP9346_A1_YP_009609871.1','PHKP9778_C1_YP_009592223.1','PH-KP9072-1_hypothesis_protein_1000_1458','PH-KP9332_A1_tr|A0A286MN37|A0A286MN37_9CAUD','PHKP9822_B1_YP_007348840.1','PH-KP6017_AUX83696.1','PHKP4049_hypothesis_protein_39750_39962','PHKP9283_B1_AZF88843.1','PHKP4049_AZF88843.1']
for i in list1:
    for key in sorted(dict1.keys()):
        if re.match(key,i):
            print(i,key,dict1[key],sep=',')