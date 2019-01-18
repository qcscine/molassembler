# TODO
# - textRepresent doesn't work quite as desired, we want 1s2 for He, not [He]
# - the valence electron count isn't right, it has to be the respective counts
#   in the current shell, meaning from the last noble gas upwards, not just non-full shells
import re

elementNameMatcher = re.compile(r"\s(?P<elementName>[A-Z][A-z]?)\s")
elementNameList = []

# read a list of element names from a file
with open("elementList.txt", "r") as inputFile:
	fullText = inputFile.read()
	for match in elementNameMatcher.finditer(fullText):
		elementNameList.append(match.group("elementName"))

vdwRadiiMatcher = re.compile(r"::(?P<elementName>[A-Z][a-z]?)\,\s(?P<number>[0-9](\.[0-9]+)?)")
vdwRadiiMap = {}

# read vdwRadii info from file
with open("vdwradii.txt", "r") as inputFile:
	fullText = inputFile.read()
	for match in vdwRadiiMatcher.finditer(fullText):
		vdwRadiiMap[str(match.group("elementName"))] = float(match.group("number"))

# comparison function for (n, l) tuples for later sorting
def nlOccupationTupleCmp(a, b):
	if a[0] < b[0]:
		return -1
	elif a[0] > b[0]:
		return 1
	else:
		if a[1] < b[1]:
			return -1
		elif a[1] > b[1]:
			return 1
		else:
			return 0

def nlOccupationTupleFull(a):
	return a[2] == 2 * (2 * a[1] + 1)

def madelung(Z):
	# generate all required orbitals for this Z
	toStore = Z
	nlOccupationTuples = []
	
	nlSum = 1
	while toStore > 0:
		# find "maximal" n, l pair to start with
		n = nlSum
		l = 0
		while l + 1 < n - 1:
			l += 1
			n -= 1

		while l >= 0:
			localStorable = 2 * (
				2 * l + 1
			)
			nlOccupationTuples.append((
				n,
				l,
				min(localStorable, toStore)
			))
			toStore -= localStorable
			if toStore <= 0: # to avoid extra 0-occupied orbitals at the end
				break
			n += 1
			l -= 1

		nlSum += 1

	return nlOccupationTuples
	#return sorted(
	#	nlOccupationTuples,
	#	cmp=nlOccupationTupleCmp
	#)


def textRepresentMadelung(madelungTupleList):
	lLetter = ["s", "p", "d", "f"]
	printstr = ""

	# find highest fully occupied (n, n-1) tuple
	highestFullyOccupiedIndex = -1
	currentShell = 1
	for c, nlOccupationTuple in enumerate(madelungTupleList):
		if nlOccupationTuple[0] == currentShell + 1 and nlOccupationTupleFull(madelungTupleList[c - 1]):
			highestFullyOccupiedIndex = c - 1
			currentShell += 1

	if highestFullyOccupiedIndex == -1: # no noble gas reached that we can reduce the representation to
		for nlOccupationTuple in madelungTupleList:
			printstr += str(nlOccupationTuple[0]) + lLetter[nlOccupationTuple[1]] + str(nlOccupationTuple[2]) + " "
	else:
		# count all electrons up to this point
		electronCountToHighest = 0
		for i in xrange(highestFullyOccupiedIndex + 1):
			electronCountToHighest += madelungTupleList[i][2]

		# Which atom is that?
		printstr += "[" + elementNameList[electronCountToHighest - 1] + "] "

		for nlOccupationTuple in madelungTupleList[highestFullyOccupiedIndex + 1:]:
			printstr += str(nlOccupationTuple[0]) + lLetter[nlOccupationTuple[1]] + str(nlOccupationTuple[2]) + " "

	return printstr

def reduceToValenceElectrons(madelungTupleList):
	# find highest fully occupied (n, n-1) tuple
	highestFullyOccupiedIndex = -1
	currentShell = 1
	for c, nlOccupationTuple in enumerate(madelungTupleList):
		if nlOccupationTuple[0] == currentShell + 1 and nlOccupationTupleFull(madelungTupleList[c - 1]):
			highestFullyOccupiedIndex = c - 1
			currentShell += 1

	maxL = 4
	counts = []
	for l in xrange(maxL):
		valenceCount = 0
		for nlOccupationTuple in madelungTupleList[highestFullyOccupiedIndex + 1:]:
			if nlOccupationTuple[1] == l:
				valenceCount += nlOccupationTuple[2]

		counts.append(valenceCount)
	
	return counts

if False:
	for c, elementName in enumerate(elementNameList):
		Z = c + 1
		madel = madelung(Z)
		valenceList = reduceToValenceElectrons(madel)
		printstr = elementName + " ("+ str(Z) + "): " + textRepresentMadelung(madel) 
		printstr += " -> " + str(valenceList) + " => " + str(sum(valenceList))

		print printstr

# make a cpp map containing the valence electrons of every shell
for c, elementName in enumerate(elementNameList):
	Z = c + 1
	madel = madelung(Z)
	valenceList = reduceToValenceElectrons(madel)
	printstr = "  {Utils::ElementType::"+elementName
	printstr += ", {" + str(vdwRadiiMap[elementName]) + ", "
	printstr += ", ".join([str(x) for x in valenceList])+"}},"

	print printstr
