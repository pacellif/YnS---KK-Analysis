import ROOT

#	OPEN ROOTUPLE

with open('Y2SPhiRun2List.txt') as f:
	allFiles = f.readlines()

for i in range(len(allFiles)):			
	allFiles[i] = allFiles[i].replace("\n", "")

def test (a = 0): 
	return allFiles if a == 0 else allFiles[:a]
