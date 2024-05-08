import csv

spreads = []
tnames = [""]

newrows = [[""]]
with open("131102 Data Spreads Set0 Centro Quad Simplified - Coord.csv") as csvfile:
    b_reader = csv.reader(csvfile)

    i = 0
    for row in b_reader:
        if i == 2:
            spread_count = 0
            row_pos = 0
            for str in row:
                if row_pos >= 2 and (row_pos - 2) % 3 == 0:                   
                    newrows.append([str]) # add a row with this spread name
                row_pos += 1
        elif i>=5:
            if row[0] != "":
                newrows[0].append(row[0] + "_x")
                newrows[0].append(row[0] + "_y")
                newrows[0].append(row[0] + "_z")
                newrows[0].append(row[0] + "_z_alt")
                newrows[0].append("")
                row_pos = 0
                for str in row:
                    if row_pos >=2 and (row_pos-2) % 3 == 0:
                        spread = int((row_pos-2)/3)
                        newrows[spread + 1].append(str)
                    elif row_pos >=3 and (row_pos-3) % 3 == 0:
                        spread = int((row_pos-3)/3)
                        newrows[spread + 1].extend([str,"","",""])
                    row_pos += 1
        i += 1

with open("b_5_oct_13.csv","w") as csvfile:
    twriter = csv.writer(csvfile)
    for row in newrows:
        twriter.writerow(row)

        



