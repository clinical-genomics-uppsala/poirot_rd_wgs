#!/bin/python3.6
import sys
import xlsxwriter
from datetime import date
import subprocess
import yaml
import gzip

# Define sys.argvs
configfile = snakemake.input[0]
duplicationFile = snakemake.input[1]
mosdepth = snakemake.input[3]  # Summary file
output = snakemake.output[0]


''' Create excel and overview tab '''
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')

# Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

with open(configfile, 'r') as file:
    config_list = yaml.load(file, Loader=yaml.FullLoader)

#runID = config_list['seqID']['sequencerun']  # sys.argv[5]
minCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[0])
medCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[1])
maxCov = int(config_list['create_cov_excel']['covLimits'].split(' ')[2])
bedfile = config_list["reference"]["coverage_bed"]
genepanels = config_list["reference"]["genepanels"]

# Define sample based on annotated vcf
sample = mosdepth.split("_")[1].split("/")[1]
today = date.today()
emptyList = ['', '', '', '', '', '']
row = 8
col = 0


''' Coverage and threshold sheet '''
covRegionsFile = mosdepth.replace(".mosdepth.summary.txt", ".regions.bed.gz")
covThresFile = mosdepth.replace('.mosdepth.summary.txt', '.thresholds.bed.gz')
tableLinesCov = []
tableLinesTre = []
totalMinBreadth = 0
totalMedBreadth = 0
totalMaxBreadth = 0
totalLength = 0
bedfile = []
covtable = []
i = 0

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
        min = int(line[4])
        med = int(line[5])
        max = int(line[6])
        totalLength += length
        totalMinBreadth += int(line[4])
        totalMedBreadth += int(line[5])
        totalMaxBreadth += int(line[6])
        row = [[line[3]],line[0:3],min,med,max,length]
        tableLinesTre.append(row)

while i < len(tableLinesCov):
    length = 0
    cov=min=med=max=0
    while (i < (len(tableLinesCov)-1)) and (tableLinesCov[i][0] == tableLinesCov[i+1][0]):
        length = length + float(tableLinesCov[i][4])
        cov = cov + (float(tableLinesCov[i][3])*float(tableLinesCov[i][4]))
        min = min + float(tableLinesTre[i][2])
        med = med + float(tableLinesTre[i][3])
        max = max + float(tableLinesTre[i][4])
        i = i + 1
    if (i == (len(tableLinesCov) - 1)) or (tableLinesCov[i][0] != tableLinesCov[i+1][0]):
        gene = tableLinesCov[i][0]
        transcript = tableLinesCov[i][1]
        length = length + float(tableLinesCov[i][4])
        cov = cov + (float(tableLinesCov[i][3])*float(tableLinesCov[i][4]))
        min = min + float(tableLinesTre[i][2])
        med = med + float(tableLinesTre[i][3])
        max = max + float(tableLinesTre[i][4])
        minBredth = min/length
        medBreadth = med/length
        maxBreadth = max/length
        covRow = [gene, transcript, round((cov/length), 2), round((minBredth*100), 1),
            round((medBreadth*100), 1), round((maxBreadth*100), 1)]
        covtable.append(covRow)
        i = i + 1


''' Low coverage sheet '''
row = 6  # 0 index
lowCovLines = []
lowRegLines = []
mosdepthPerBase = mosdepth.replace(".mosdepth.summary.txt", ".mosdepth.lowCov.regions.txt")

with open(mosdepthPerBase, 'r') as file:
    for lline in file:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(minCov):
            lowCovLines.append(line)

for line in lowCovLines:
    for bedline in bedfile:
        if line[0] == bedline[0] and int(line[1]) >= int(bedline[1]) and int(line[2]) <= int(bedline[2]):
            lowregion=[bedline[3],line[0],line[1],line[2],line[3]]
            lowRegLines.append(lowregion)
        elif line[0] == bedline[0] and int(line[1]) < int(bedline[1]) and int(line[2]) >= int(bedline[1]):
            if int(line[2]) >= int(bedline[2]):
                lowregion=[bedline[3],line[0],bedline[1],bedline[2],line[3]]
                lowRegLines.append(lowregion)
                break
            else:
                lowregion=[bedline[3],line[0],bedline[1],line[2],line[3]]
                lowRegLines.append(lowregion)
        elif line[0] == bedline[0] and int(line[1]) < int(bedline[2]) and int(line[2]) > int(bedline[2]):
            lowregion=[bedline[3],line[0],line[1],bedline[2],line[3]]
            lowRegLines.append(lowregion)
            break


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

# Add avg. cov and clonality
cmdAvgCov = "grep total_region "+mosdepth+" | awk '{print $4}'"
avgCov = subprocess.run(cmdAvgCov, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

cmdDupl = 'grep -A1 PERCENT '+duplicationFile+' | tail -1 | cut -f9'
duplicateLevel = subprocess.run(cmdDupl, stdout=subprocess.PIPE, shell='TRUE').stdout.decode('utf-8').strip()

worksheetOver.write_row(8, 0, ['DNAnr', 'Avg. coverage (x)', 'Duplicationlevel ()%)',
                                str(minCov)+'x (%)', str(medCov)+'x (%)', str(maxCov)+'x (%)'], tableHeadFormat)
worksheetOver.write_row(9, 0, [sample, avgCov, str(round(float(duplicateLevel)*100, 2)),
                                str(round(float(totalMinBreadth/totalLength)*100, 1)),
                                str(round(float(totalMedBreadth/totalLength)*100, 1)),
                                str(round(float(totalMaxBreadth/totalLength)*100, 1))])

worksheetOver.write(11, 0, 'Number of regions not covered by at least '+str(minCov)+'x: ', tableHeadFormat)  # From Cov sheet
worksheetOver.write(12, 0, str(len(lowRegLines)))  # From Cov sheet

worksheetOver.write(14, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(15, 0, "internal:'Low Coverage'!A1", string='Positions with coverage lower than '+str(minCov)+'x')
worksheetOver.write_url(16, 0, "internal: 'Coverage'!A1", string='Average coverage of all regions in bed')
#worksheetOver.write_url(10, 0, "internal:'Version'!A1", string='Version Log')


''' Genepanel sheets '''
panels = []
with open(genepanels, 'rt') as file:
    for lline in file:
        panels.append(lline.rstrip())

number = 0
for line in panels:
    panel = genepanels.replace("genepanels", line)
    worksheetpanel = workbook.add_worksheet(line)
    worksheetpanel.write('A1', 'Coverage analysis per gene for '+line+' gene panel', headingFormat)
    worksheetpanel.write_row('A2', emptyList, lineFormat)
    worksheetpanel.write('A3', 'Sample: '+str(sample))
    worksheetpanel.write('A5', 'Averge coverage and coverage breadth of genes in '+line+' gene panel.')
    number = number + 1
    genes = []
    lows = []
    with open(panel, 'rt', encoding='utf-8') as file:
        for lline in file:
            i = 0
            k = 0
            while i < len(covtable):
                if lline.rstrip() == covtable[i][0]:
                    genes.append(covtable[i])
                    break
                else:
                    i = i + 1
            while k < len(lowRegLines):
                lowgene = lowRegLines[k][0].strip().split("_")
                if lline.rstrip() == lowgene[0]:
                    lows.append(lowRegLines[k])
                    k = k + 1
                else:
                    k = k + 1
        tableArea = 'A6:F'+str(len(genes)+6)  # rows of full list
        headerListDict = [{'header': 'Gene'}, {'header': 'Transcript'}, {'header': 'Avg coverage'},
                            {'header': str(minCov)+'x'}, {'header': str(medCov)+'x'}, {'header': str(maxCov)+'x'}]
        worksheetpanel.add_table(tableArea, {'data': genes, 'columns': headerListDict, 'style': 'Table Style Light 1'})
        if lows != []:
            worksheetpanel.write('A'+str(len(genes)+9), 'Regions of exons that are covered below '+str(minCov)+'x'+'.')
            tableArea2 = 'A'+str(len(genes)+10)+':E'+str(len(lows)+len(genes)+10)  # rows of full list
            headerListDict2 = [{'header': 'Gene name_transcript_exon'}, {'header': 'Chr'}, {'header': 'Start'},
                                {'header': 'Stop'}, {'header': 'Coverage (x)'}]
            worksheetpanel.add_table(tableArea2, {'data': lows, 'columns': headerListDict2, 'style': 'Table Style Light 1'})
        else:
            worksheetpanel.write('A'+str(len(genes)+9), 'No regions of exons are covered under '+str(minCov)+'x'+'.')
        worksheetOver.write_url((16+number), 0, "internal: '"+line+"'!A1", string=line)


''' Coverage sheet '''
worksheetCov = workbook.add_worksheet('Coverage')
worksheetCov.write('A1', 'Average coverage and coverage breadth per gene', headingFormat)
worksheetCov.write_row('A2', emptyList, lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A5', 'Averge coverage and coverage breadth of each gene in exon-bedfile')
tableArea = 'A6:F'+str(len(covtable)+6)  # rows of full list
headerListDict = [{'header': 'Gene'}, {'header': 'Transcript'}, {'header': 'Avg coverage'},
                    {'header': str(minCov)+'x'}, {'header': str(medCov)+'x'}, {'header': str(maxCov)+'x'}]
worksheetCov.add_table(tableArea, {'data': covtable, 'columns': headerListDict, 'style': 'Table Style Light 1'})


''' Low Coverage sheet '''
worksheetLowCov = workbook.add_worksheet('Low Coverage')
worksheetLowCov.write('A1', 'Mosdepth coverage analysis', headingFormat)
worksheetLowCov.write_row('A2', emptyList, lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
worksheetLowCov.write('A5', 'Gene regions with coverage lower than '+str(minCov)+'x.')
tableArea = 'A6:E'+str(len(lowRegLines)+6)  # rows of full list
headerListDict = [{'header': 'Region Name'}, {'header': 'Chr'}, {'header': 'Start'},
                    {'header': 'Stop'}, {'header': 'Mean Coverage'}]
worksheetLowCov.add_table(tableArea, {'data': lowRegLines, 'columns': headerListDict, 'style': 'Table Style Light 1'})


workbook.close()
