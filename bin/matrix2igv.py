#!/usr/bin/env python
import csv
import argparse
import gzip

sampleEndMR = {}
sample_names = []
posSample = {}

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="specify the input matrix file", required=True)
    parser.add_argument("-o", "--output", help="specify the output file", required=True)
    parser.add_argument("-m", "--meth", help = "specify the input methylated region file", required=True)
    args = parser.parse_args()

    createCSV(args, parser)

    #command = "igvtools toTDF --fileType igv -t . " + args.output + " " + args.output.replace(".igv", ".tdf")
    #print command
    #os.system(command)


def createCSV(args, parser):

    newcsvfile = open(args.output, 'w')
    write_in_file = csv.writer(newcsvfile, delimiter='\t')
    newcsvfile.write("#track viewLimits=0:100 coords=0\n")

    if args.input.endswith(".gz"):
        matrixFile = gzip.open(args.input, 'r')
    else:
        matrixFile = open(args.input, 'r')

    matrixFile = csv.reader(matrixFile, delimiter='\t')

    with open(args.meth, 'r' ) as mrFile:
        check_header = True

        #DefaultLine to start procedure in iterThroughMR
        actualMRLine = "defaultLine"

        #Dealing with the samples
        for lineMatrix in matrixFile:

            if len(lineMatrix) <= 1:
                continue

            #create the header line, posSample and sampleEndMR
            # posSample = hashmap to get the belonging Samplename of the value in row
            # sampleEndMR = hashmap to save every Sample with an MR at the actual write
            # position, and to get their next methylation endingpoint
            if lineMatrix[0].startswith('#'):
                file_meth_rate_row = createHeader(0, lineMatrix)
                check_header = False

            #Dealing with the samples
            else:
                # goUntil: saves the actual position
                if check_header:
                    parser.error("matrix file needs a header like: '#chr    pos    class    "
                                 "strand    sample1    sample2'")

                goUntil = {'chr': lineMatrix[0], 'pos': lineMatrix[1]}
                actualMRLine = iterThroughMR(mrFile, goUntil, actualMRLine, parser)

                file_meth_rate_row = [lineMatrix[0], int(lineMatrix[1])-1, lineMatrix[1],
                                      lineMatrix[2]]

                for i, sample in enumerate(lineMatrix[4:]):
                    #get the Sample name the value belongs to
                    #actualSamp = posSample[lineMatrix.index(sample)]
                    actualSamp = posSample[i+4]
                    # check if sample belongs to an Sample with an MR right now
                    if actualSamp in sampleEndMR:
                        if int(sampleEndMR[actualSamp]) >= int(goUntil['pos']):
                            file_meth_rate_row = getSampleInMR(sample, file_meth_rate_row)
                        #If end of the MR is reached, remove sample of the hashmap
                        else:
                            del sampleEndMR[actualSamp]
                            file_meth_rate_row = getSampleNormal(sample, file_meth_rate_row)
                    else:
                        file_meth_rate_row = getSampleNormal(sample, file_meth_rate_row)

            write_in_file.writerow(file_meth_rate_row)

    newcsvfile.close()


def getSampleInMR(sample, row):
    if not sample == '0':
        items = sample.split("/")
        if(items[1] ==0):
            row.append(0)
        else:
            row.append(int(float(items[1])))
        row.append(int(items[2])+int(items[3]))
    else:
        row.append(0)
        row.append(0)
    return row

def getSampleNormal(sample, row):
    if not sample == '0':
        items = sample.split("/")
        if(float(items[1]) == 0.0):
            row.append(0)
        else:
            row.append(int(float(items[1])*-1))
        row.append(int(items[2])+int(items[3]))
    else:
        row.append(0)
        row.append(0)
    return row


def createHeader(counter, matrixHeader):
    for i, tab in enumerate(matrixHeader):
        if i>3:
            posSample.update({4+counter: tab})
            posSample.update({5+counter: tab})
            counter += 2
            sample_names.append(tab.replace("EPI_", "") + '_rate')
            sample_names.append(tab.replace("EPI_", "") + '_cov')
    file_meth_rate_row = ['Chromosome', 'Start', 'End', 'Feature']
    file_meth_rate_row.extend(sample_names)
    return(file_meth_rate_row)


def iterThroughMR(mrFile, goUntil, actualMRLine, parser):

    mrTabFile = csv.reader(mrFile, delimiter ='\t')
    try:
        if actualMRLine == "defaultLine":
           actualMRLine = next(mrFile).split('\t')

        if actualMRLine is None:
            return actualMRLine

        if int(actualMRLine[1]) >= int(goUntil['pos']) or actualMRLine[0] != goUntil['chr']:
            return actualMRLine
        else:
            sampleEndMR.update({actualMRLine[7]: actualMRLine[2]})
            for i, line in enumerate(mrTabFile):
                if(len(line) < 1):
                    parser.error("Metyhlated regions file should not contain empty lines")
                if (int(line[1])<=int(goUntil['pos'])):
                    sampleEndMR.update({line[7]: line[2]})
                else:
                    return line
    except IndexError:
        parser.error("Methylated region file has incorrect format. No empty lines "
                     "allowed. Needed format: 'Chromosome, start, end, ... , SampleName")

if __name__ == '__main__':
    main()
