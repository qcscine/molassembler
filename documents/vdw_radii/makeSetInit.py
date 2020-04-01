import csv

with open("concat.csv", "r") as csvfile:
	csvreader = csv.reader(csvfile)

	for line in csvreader:
		print("    {Utils::ElementType::"+line[2]+", "+line[3]+"},")

