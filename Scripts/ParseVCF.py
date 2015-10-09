"""This is a custom scripts used to parse a vcf file and pull out relavant information as well as the allele state

Written by Jean-Michel Michno"""

seqfile = open("FNOnethruElevenv2Homo.vcf", "r")


# make a list of dictionaries to store all of our
Chromosome = list()
Position = list()
ID = list()
Ref_Base = list()
Alt_Base = list()
Ref = list()
Quality = list()
FilterInfo = list()
Format = list()
FN01 = list()
FN02 = list()
FN03 = list()
FN04 = list()
FN05 = list()
FN06 = list()
FN07 = list()
FN08 = list()
FN09 = list()
FN10 = list()
FN11 = list()

# make a list for the #coments and the actual data
identifierlist = list()
sequenceinfolist = list()

for line in seqfile:
    line = line.rstrip()
    if line.startswith("#"):
        identifierlist.append(line)
    else:
        sequenceinfolist.append(line)

# pull out the header info
header = identifierlist[-1]

# parse our the file
for line in sequenceinfolist:
    columns = line.split("\t")
    Chromosome.append(columns[0])
    Position.append(columns[1])
    ID.append(columns[2])
    Ref_Base.append(columns[3])
    Alt_Base.append(columns[4])
    Ref.append(columns[5])
    Quality.append(columns[6])
    FilterInfo.append(columns[7])
    Format.append(columns[8])
    FN01.append(columns[9][0:3])
    FN02.append(columns[10][0:3])
    FN03.append(columns[11][0:3])
    FN04.append(columns[12][0:3])
    FN05.append(columns[13][0:3])
    FN06.append(columns[14][0:3])
    FN07.append(columns[15][0:3])
    FN08.append(columns[16][0:3])
    FN09.append(columns[17][0:3])
    FN10.append(columns[18][0:3])
    FN11.append(columns[19][0:3])

output = open('FNOnethruElevenv2Homo.vcf', 'w')

output.write('Chromosome'+', '+'Position'+', '+'ID'+', '+'Ref_Base'+
    ', '+'Alt_Base'+', '+'Ref'+', '+'Quality'+', '+'FilterInfo'+', '+'Format'+', '+'FN01'+', '+'FN02'+', '+'FN03'+', '+'FN04'+', '+'FN05'+', '+'FN06'+', '+'FN07'+', '+'FN08'+', '+'FN09'+', '+'FN10'+', '+'FN11'+'\n')
for a, b, c, d, e, f, g, h, i, j, k, l, m, n, l, o, p, q, r, s in zip(
    Chromosome, Position, ID, Ref_Base, Alt_Base, Ref,
    Quality, FilterInfo, Format, FN01, FN02,
    FN03, FN04, FN05, FN06, FN07, FN08, FN09,
        FN10, FN11):

    final = (a, b, c, d, e, f, g, h, i, j, k, l, m, n, l, o, p, q, r, s)
    output.write(','.join(final))
    output.write('\n')

output.close()
seqfile.close()
