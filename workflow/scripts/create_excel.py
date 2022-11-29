#!/bin/python3.6
import sys
import xlsxwriter
from datetime import date
import subprocess
import yaml
import gzip

# Define sys.argvs
mosdepth = snakemake.input[0]  # Summary file
output = snakemake.output[0]
configfile = snakemake.input[2]


with open(configfile, 'r') as file:
    config_list = yaml.load(file, Loader=yaml.FullLoader)

#runID = config_list['seqID']['sequencerun']  # sys.argv[5]
minCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[0])
medCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[1])
maxCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[2])
bedfile = config_list["reference"]["coverage_bed"]


''' Create execl file and sheets. '''
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')
worksheetLowCov = workbook.add_worksheet('Low Coverage')
worksheetCov = workbook.add_worksheet('Coverage')
worksheetThers = workbook.add_worksheet('Thresholds')
# worksheetVersions = workbook.add_worksheet('Version')
# Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

# Define sample based on annotated vcf
sample = mosdepth.split("_")[1].split("/")[1]
today = date.today()
emptyList = ['', '', '', '', '', '']

# ''' Prog Version sheet '''
# worksheetVersions.write('A1', 'Version Log', headingFormat)
# worksheetVersions.write_row(1, 0, emptyList, lineFormat)
# worksheetVersions.write('A3', 'Sample: '+str(sample))
#
# worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)
# containers = [clist for clist in config_list['singularitys'].items()]
row = 8
col = 0
# for containerTuple in containers:
#     container = list(containerTuple)
#     worksheetVersions.write_row('A'+str(row), container)
#     row += 1


''' Coverage '''
covRegionsFile = mosdepth.replace(".mosdepth.summary.txt", ".regions.bed.gz")

worksheetCov.write('A1', 'Average Coverage per exon', headingFormat)
worksheetCov.write_row('A2', emptyList, lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A4', 'Averge coverage of each region in exon-bedfile')

tableLines = []
bedfile = []
with gzip.open(covRegionsFile, 'rt') as regionsfile:
    for lline in regionsfile:
        line = lline.strip().split('\t')
        gene = line[3].split("_")[0]
        exon = line[3].split("_")[3]
        transcript = "NM_"+line[3].split("_")[2]
        covRow = [line[0], line[1], line[2], gene, exon, transcript, line[4]]
        tableLines.append(covRow)
        bedfile.append(line[0:4])

tableArea = 'A6:G'+str(len(tableLines)+6)  # rows of full list
headerListDict = [{'header': 'Chr'}, {'header': 'Start'}, {'header': 'Stop'},
                  {'header': 'Gene'}, {'header': 'Exon'}, {'header': 'Transcript'}, {'header': 'Avg coverage'}, ]
worksheetCov.add_table(tableArea, {'data': tableLines, 'columns': headerListDict, 'style': 'Table Style Light 1'})


''' Threshold sheet '''
worksheetThers.set_column(1, 3, 10)
worksheetThers.set_column(1, 4, 10)
worksheetThers.write('A1', 'Coverage breadth per exon', headingFormat)
worksheetThers.write_row('A2', emptyList, lineFormat)
worksheetThers.write('A3', 'Sample: '+str(sample))
worksheetThers.write('A4', 'Coverage breadth of each region in exon-bedfile')
covThresFile = mosdepth.replace('.mosdepth.summary.txt', '.thresholds.bed.gz')

tableLines = []
totalMinBreadth = 0
totalMedBreadth = 0
totalMaxBreadth = 0
totalLength = 0
with gzip.open(covThresFile, 'rt') as thresfile:
    next(thresfile)
    for lline in thresfile:
        line = lline.strip().split("\t")
        length = int(line[2])-int(line[1])
        totalLength += length
        totalMinBreadth += int(line[4])
        totalMedBreadth += int(line[5])
        totalMaxBreadth += int(line[6])
        minBredth = round(int(line[4])/length, 4)
        medBreadth = round(int(line[5])/length, 4)
        maxBreadth = round(int(line[6])/length, 4)
        tableLines.append([line[3]]+line[0:3]+[str(minBredth)]+[str(medBreadth)]+[str(maxBreadth)])

tableArea = 'A6:G'+str(len(tableLines)+6)  # rows of full list
headerListDict = [{'header': 'Regions'}, {'header': 'Chr'}, {'header': 'Start'}, {'header': 'Stop'}, {'header': str(minCov)+'x'},
                  {'header': str(medCov)+'x'}, {'header': str(maxCov)+'x'}]
worksheetThers.add_table(tableArea, {'data': tableLines, 'columns': headerListDict, 'style': 'Table Style Light 1'})

''' Low Coverage sheet '''
worksheetLowCov.set_column(1, 3, 10)
worksheetLowCov.set_column(1, 4, 10)
# Heading in sheet
worksheetLowCov.write('A1', 'Mosdepth coverage analysis', headingFormat)
worksheetLowCov.write_row('A2', emptyList, lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetLowCov.write('A4', description)
covHeadings = ['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage']
worksheetLowCov.write_row('A6', covHeadings, tableHeadFormat)  # 1 index
row = 6  # 0 index

lowCov = []
mosdepthPerBase = mosdepth.replace(".mosdepth.summary.txt", ".mosdepth.lowCov.regions.txt")

lowCovLines = []
with open(mosdepthPerBase, 'r') as file:
    for lline in file:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(minCov):
            lowCovLines.append(line)

# Sort based on coverage
lowCovLines.sort(key=lambda x: x[3])
for line in lowCovLines:
    for bedline in bedfile:
        if line[0] == bedline[0] and int(line[1]) >= int(bedline[1]) and int(line[2]) <= int(bedline[2]):
            worksheetLowCov.write_row(row, col, [bedline[3]]+line)
            row += 1
            break

# Number of low cov regions for the Overview sheet.
lowRegions = row - 6


''' Overview sheet (1) '''
worksheetOver.write(0, 0, sample, headingFormat)
#worksheetOver.write(1, 0, "RunID: "+runID)
worksheetOver.write(2, 0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(3, 0, emptyList, lineFormat)

worksheetOver.write(4, 0, "Created by: ")
worksheetOver.write(4, 4, "Valid from: ")
worksheetOver.write(5, 0, "Signed by: ")
worksheetOver.write(5, 4, "Document nr: ")
worksheetOver.write_row(6, 0, emptyList, lineFormat)

worksheetOver.write(7, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(8, 0, "internal:'Low Coverage'!A1", string='Positions with coverage lower than '+str(minCov)+'x')
worksheetOver.write_url(9, 0, "internal: 'Coverage'!A1", string='Average coverage of all regions in bed')
#worksheetOver.write_url(10, 0, "internal:'Version'!A1", string='Version Log')
worksheetOver.write_row(11, 0, emptyList, lineFormat)
worksheetOver.write_row(10, 0, emptyList, lineFormat)


# # Add avg. cov and clonality
# cmdAvgCov = "grep total_region "+mosdepth+" | awk '{print $4}'"
# avgCov = subprocess.run(cmdAvgCov, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
#
# duplicationFile = mosdepth.replace(".mosdepth.summary.txt", "_DuplicationMetrics.txt")
# cmdDupl = 'grep -A1 PERCENT '+duplicationFile+' | tail -1 | cut -f9'
# duplicateLevel = subprocess.run(cmdDupl, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
#
#
# worksheetOver.write_row(13, 0, ['RunID', 'DNAnr', 'Avg. coverage [x]', 'Duplicationlevel [%]',
#                                 str(minCov)+'x', str(medCov)+'x', str(maxCov)+'x'], tableHeadFormat)
# worksheetOver.write_row(14, 0, [runID, sample, avgCov, str(round(float(duplicateLevel)*100, 2)),
#                                 str(round(float(totalMinBreadth/totalLength), 4)),
#                                 str(round(float(totalMedBreadth/totalLength), 4)),
#                                 str(round(float(totalMaxBreadth/totalLength), 4))])
#
#
# worksheetOver.write(16, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ')  # From Cov sheet
# worksheetOver.write(17, 0, str(lowRegions))  # From Cov sheet

worksheetOver.write(11, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ')  # From Cov sheet
worksheetOver.write(12, 0, str(lowRegions))  # From Cov sheet


workbook.close()







#!/bin/python3.6
import sys
import xlsxwriter
from datetime import date
import subprocess
import yaml
import gzip

# Define sys.argvs
mosdepth="qc/mosdepth_bed/29893_N.mosdepth.summary.txt"
output = "qc/mosdepth_bed/29893_N.coverage.xlsx"
configfile = "config/config.yaml"


with open(configfile, 'r') as file:
    config_list = yaml.load(file, Loader=yaml.FullLoader)

#runID = config_list['seqID']['sequencerun']  # sys.argv[5]
minCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[0])
medCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[1])
maxCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[2])
bedfile = config_list["reference"]["coverage_bed"]


''' Create execl file and sheets. '''
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')
worksheetLowCov = workbook.add_worksheet('Low Coverage')
worksheetCov = workbook.add_worksheet('Coverage')
worksheetThers = workbook.add_worksheet('Thresholds')
# worksheetVersions = workbook.add_worksheet('Version')
# Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

# Define sample based on annotated vcf
sample = mosdepth.split("_")[1].split("/")[1]
today = date.today()
emptyList = ['', '', '', '', '', '']

# ''' Prog Version sheet '''
# worksheetVersions.write('A1', 'Version Log', headingFormat)
# worksheetVersions.write_row(1, 0, emptyList, lineFormat)
# worksheetVersions.write('A3', 'Sample: '+str(sample))
#
# worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)
# containers = [clist for clist in config_list['singularitys'].items()]
 row = 8
 col = 0
# for containerTuple in containers:
#     container = list(containerTuple)
#     worksheetVersions.write_row('A'+str(row), container)
#     row += 1


''' Coverage and threshold sheet '''
covRegionsFile = mosdepth.replace(".mosdepth.summary.txt", ".regions.bed.gz")
covThresFile = mosdepth.replace('.mosdepth.summary.txt', '.thresholds.bed.gz')
worksheetCov.write('A1', 'Average coverage and coverage breadth per gene', headingFormat)
worksheetCov.write_row('A2', emptyList, lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A4', 'Averge coverage and coverage breadth of each gene in exon-bedfile')

tableLinesCov = []
tableLinesTre = []
totalMinBreadth = 0
totalMedBreadth = 0
totalMaxBreadth = 0
totalLength = 0
bedfile = []

with gzip.open(covRegionsFile, 'rt') as regionsfile:
    for lline in regionsfile:
        line = lline.strip().split('\t')
        length = int(line[2])-int(line[1])
        gene = line[3].split("_")[0]
        exon = line[3].split("_")[3]
        transcript = "NM_"+line[3].split("_")[2]
        covRow = [gene, transcript, exon, line[4], length]
        tableLinesCov.append(covRow)
        bedfile.append(line[0:4])

with gzip.open(covThresFile, 'rt') as thresfile:
    next(thresfile)
    for lline in thresfile:
        line = lline.strip().split("\t")
        length = int(line[2])-int(line[1])
        totalLength += length
        min = int(line[4])
        med = int(line[4])
        max = int(line[4])
        totalMinBreadth += int(line[4])
        totalMedBreadth += int(line[5])
        totalMaxBreadth += int(line[6])
        row = [[line[3]],line[0:3],Min,Med,Max,length]
        tableLinesTre.append(row)

        minBredth = int(line[4])/length
        medBreadth = int(line[5])/length
        maxBreadth = int(line[6])/length
covtable = []
i = 0
while i <= (len(tableLinesCov)-1):
	length = 0
	cov=min=med=max=0
	while (i < (len(tableLinesCov)-1)) and (tableLinesCov[i][0] == tableLinesCov[i+1][0]):
		length = length + float(tableLinesCov[i][4])
		cov = cov + (float(tableLinesCov[i][3])*float(tableLinesCov[i][4]))
		min = min + float(tableLinesTre[i][4])
		med = med + float(tableLinesTre[i][5])
		max = max + float(tableLinesTre[i][6])
		i = i + 1
	if (i == (len(tableLinesCov) - 1)) or (tableLinesCov[i][0] != tableLinesCov[i+1][0]):
		gene = tableLinesCov[i][0]
		transcript = tableLinesCov[i][1]
		length = length + float(tableLinesCov[i][4])
		cov = cov + (float(tableLinesCov[i][3])*float(tableLinesCov[i][4]))
		min = min + float(tableLinesTre[i][4])
		med = med + float(tableLinesTre[i][5])
		max = max + float(tableLinesTre[i][6])
		covRow = [gene, transcript, round((cov/length), 2), round((min/exon), 2), round((med/exon), 2), round((max/exon), 2)]
		covtable.append(covRow)
		i = i + 1

tableArea = 'A6:F'+str(len(covtable)+6)  # rows of full list
headerListDict = [{'header': 'Gene'}, {'header': 'Transcript'}, {'header': 'Avg coverage'},
	{'header': str(minCov)+'x'}, {'header': str(medCov)+'x'}, {'header': str(maxCov)+'x'}]
worksheetCov.add_table(tableArea, {'data': covtable, 'columns': headerListDict, 'style': 'Table Style Light 1'})


''' Low Coverage sheet '''
worksheetLowCov.set_column(1, 3, 10)
worksheetLowCov.set_column(1, 4, 10)
# Heading in sheet
worksheetLowCov.write('A1', 'Mosdepth coverage analysis', headingFormat)
worksheetLowCov.write_row('A2', emptyList, lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetLowCov.write('A4', description)
covHeadings = ['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage']
worksheetLowCov.write_row('A6', covHeadings, tableHeadFormat)  # 1 index
row = 6  # 0 index

lowCov = []
mosdepthPerBase = mosdepth.replace(".mosdepth.summary.txt", ".mosdepth.lowCov.regions.txt")

lowCovLines = []
with open(mosdepthPerBase, 'r') as file:
    for lline in file:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(minCov):
            lowCovLines.append(line)

# Sort based on coverage
lowCovLines.sort(key=lambda x: x[3])
for line in lowCovLines:
    for bedline in bedfile:
        if line[0] == bedline[0] and int(line[1]) >= int(bedline[1]) and int(line[2]) <= int(bedline[2]):
            worksheetLowCov.write_row(row, col, [bedline[3]]+line)
            row += 1
            break

# Number of low cov regions for the Overview sheet.
lowRegions = row - 6


''' Overview sheet (1) '''
worksheetOver.write(0, 0, sample, headingFormat)
worksheetOver.write(1, 0, "RunID: "+runID)
worksheetOver.write(2, 0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(3, 0, emptyList, lineFormat)

worksheetOver.write(4, 0, "Created by: ")
worksheetOver.write(4, 4, "Valid from: ")
worksheetOver.write(5, 0, "Signed by: ")
worksheetOver.write(5, 4, "Document nr: ")
worksheetOver.write_row(6, 0, emptyList, lineFormat)

worksheetOver.write(7, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(8, 0, "internal:'Low Coverage'!A1", string='Positions with coverage lower than '+str(minCov)+'x')
worksheetOver.write_url(9, 0, "internal: 'Coverage'!A1", string='Average coverage of all regions in bed')
#worksheetOver.write_url(10, 0, "internal:'Version'!A1", string='Version Log')
worksheetOver.write_row(11, 0, emptyList, lineFormat)
worksheetOver.write_row(10, 0, emptyList, lineFormat)


# Add avg. cov and clonality
cmdAvgCov = "grep total_region "+mosdepth+" | awk '{print $4}'"
avgCov = subprocess.run(cmdAvgCov, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

# duplicationFile = mosdepth.replace(".mosdepth.summary.txt", "_DuplicationMetrics.txt")
# cmdDupl = 'grep -A1 PERCENT '+duplicationFile+' | tail -1 | cut -f9'
# duplicateLevel = subprocess.run(cmdDupl, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()
#
#
# worksheetOver.write_row(13, 0, ['RunID', 'DNAnr', 'Avg. coverage [x]', 'Duplicationlevel [%]',
#                                 str(minCov)+'x', str(medCov)+'x', str(maxCov)+'x'], tableHeadFormat)
# worksheetOver.write_row(14, 0, [sample, avgCov, str(round(float(duplicateLevel)*100, 2)),
#                                 str(round(float(totalMinBreadth/totalLength), 4)),
#                                 str(round(float(totalMedBreadth/totalLength), 4)),
#                                 str(round(float(totalMaxBreadth/totalLength), 4))])

worksheetOver.write_row(14, 0, [sample, avgCov,
                                str(round(float(totalMinBreadth/totalLength), 4)),
                                str(round(float(totalMedBreadth/totalLength), 4)),
                                str(round(float(totalMaxBreadth/totalLength), 4))])
#
#
# worksheetOver.write(16, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ')  # From Cov sheet
# worksheetOver.write(17, 0, str(lowRegions))  # From Cov sheet

worksheetOver.write(11, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ')  # From Cov sheet
worksheetOver.write(12, 0, str(lowRegions))  # From Cov sheet


workbook.close()
